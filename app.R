require(shiny)

options(shiny.maxRequestSize = 30*1024^2)

library(data.table)
library(DT)
library(abind)
library(scater)
library(SummarizedExperiment)
library(RColorBrewer)
library(shiny)
library(pheatmap)
library(reader)

fPem <- function(x)
{
    if(!all(is.na(x)))
    {
        x <- as.matrix(x)
        x[x<0] <- NA
        x <- cbind(x, r=rowSums(x, na.rm=FALSE)) 
        x <- rbind(x, c=colSums(x, na.rm=TRUE))	
        x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <-
            x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] /
            (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x),
                                                    which(x[nrow(x),]>0)] /
                 x[nrow(x), ncol(x)]) 
        
        x <- log10(x)
        
    } else {
        x <- NA
        print("No data avalable.")
    }
    
    x[which(!is.finite(x))] <- 0
    
    return(x)
}

shinyApp(
    
    ui = navbarPage(
        
        title = "Preferential Expression Measure Calculator",
        id="Preferential Expression Measure Calculator",
        fluid=TRUE,
        
        tabPanel("Upload Tab",
                 sidebarLayout(
                     sidebarPanel(
                         radioButtons("uploadChoice", "",
                                      c("Count File and Metadata File" = "countFile",
                                        "Summarized Experiment Object" = "seObject"
                                      )),
                         conditionalPanel(condition = "input.uploadChoice == 'seObject'",
                                          fileInput(
                                              "se",
                                              "Upload Summarized Experiment",
                                              multiple = FALSE,
                                              accept = ".RDS"
                                          )),
                         conditionalPanel(condition = "input.uploadChoice == 'countFile'",
                                          fileInput(
                                              "counts",
                                              "Upload Counts Table",
                                              multiple = FALSE,
                                              accept = ".csv"
                                          ),
                                          fileInput("md", "Upload Metadata Table",
                                                    multiple = FALSE,
                                                    accept = ".csv")),
                         actionButton('uploads', label = 'Upload'),
                     ),
                     mainPanel(
                         dataTableOutput("summary"),
                     )
                 )
        ),
        tabPanel("PEM Analysis Tab",
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput('assay','Choose assay for analysis',
                                        multiple=FALSE,
                                        choices = c(''),
                                        selected = NULL),
                         selectizeInput('condition','Choose metadata for analysis grouping',
                                        multiple=FALSE,
                                        choices = c(''),
                                        selected = NULL),
                         actionButton('analyze', label = 'Analyze PEM'),
                     ),
                     mainPanel(
                         dataTableOutput("pem"),
                         conditionalPanel(condition = "output.pem",
                                          downloadButton("download_pems", "Download")),
                     )
                 )
        ),
        tabPanel("Visualization Tab",
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput('plot_type','Choose visualization',
                                        multiple=FALSE,
                                        choices = c('scatter','bar'),
                                        selected = NULL),
                         selectizeInput('gene_name','Choose gene to visualize',
                                        multiple=FALSE,
                                        choices = c(''),
                                        selected = NULL),
                         actionButton('plot', label = 'Plot')
                     ),
                     mainPanel(
                         plotOutput('plots'),
                         conditionalPanel(condition = "output.tab",
                                          p('Percents of PEM Scores per Metadata Group')),
                         dataTableOutput("tab")
                     )
                 )
        ),
    ), 
    server = function(input, output) {
        
        reactivevalue=reactiveValues(se=NULL,
                                     data=NULL,
                                     scores=NULL
        )
        
        observeEvent(input$uploads, {
            
            if (input$uploadChoice == 'countFile') {
                
            counts_matrix <- input$counts$datapath
            
            matrix <- read.table(counts_matrix,
                             header = TRUE, row.names = 1,
                             sep = get.delim(counts_matrix,
                                             n = 10,
                                             delims = c('\t',',')))
            
            metadata <- input$md$datapath
            
            md <- read.table(metadata,
                                                 header = TRUE, row.names = 1,
                                                 sep = get.delim(metadata,
                                                                 n = 10,
                                                                 delims = c('\t',',')))
            
            se <- SummarizedExperiment(list(counts = matrix),
                                       colData = md)
            reactivevalue$se <- se
            
            }
            
            if (input$uploadChoice == 'seObject') {
            
            se <- readRDS(input$se$datapath)
            
            reactivevalue$se <- se
            
            }
            
            updateSelectizeInput(inputId = "assay",
                                 choices = assayNames((reactivevalue$se)),
                                 selected = NULL)
            
            updateSelectizeInput(inputId = "condition",
                                 choices = colnames(colData(reactivevalue$se)),
                                 selected = NULL)
            
            md <- as.data.frame(colData(reactivevalue$se))
            
            md[is.na(md)] = 0
            CN <- colnames(md)
            Type <- sapply(md, class)
            SN <- c()
            MDV <- c()
            for (i in 1:length(colnames(md))) {
                if (Type[i]=="numeric") {
                    st <- paste("+/-", round(sd(md[,i])))
                    MDV <- c(MDV,paste(round(mean(md[,i])), st))
                    SN <- c(SN,'N/A')
                }
                else if (Type[i]=="integer") {
                    st <- paste("+/-", round(sd(md[,i])))
                    MDV <- c(MDV,paste(round(mean(md[,i])), st))
                    SN <- c(SN,'N/A')
                }
                else {
                    MDV <- c(MDV,paste(unique(md[,i]), sep=", ", collapse=", "))
                    n <- c()
                    for (j in unique(md[,i])) {
                        n <- c(n,length(which(md[,i]==j)))
                    }
                    SN <- c(SN,paste(n,collapse = '/'))
                }
            }
            sum_tab <- abind(CN,Type,MDV,SN,along=2)
            col_names <- c('Column Name','Type','Mean (sd) or Distinct Values','Sample Count')
            colnames(sum_tab) <- col_names
            output$summary <- renderDataTable((as.data.table(sum_tab)))
        }
        ) 
        
        observeEvent(input$analyze, {
            
            withProgress({
                setProgress(0.5, 'Calculating...')
                count_data <- as.data.frame(reactivevalue$se@assays@data[[input$assay]])
                
                coldata <- colData(reactivevalue$se)
                
                colnames(count_data) <- coldata[,input$condition]
                
                mean_data <- as.data.frame(sapply(unique(names(count_data)),
                                                  function(col)
                                                      rowMeans(
                                                          count_data[names(
                                                              count_data) == col]) 
                )
                )
                
                rownames(mean_data) <- rownames(reactivevalue$se)
                
                reactivevalue$data <- mean_data
                
                pems <- fPem(mean_data)
                
                pems <- pems[-c(which(rownames(pems)=="c")),]
                
                pems <- pems[,-c(which(colnames(pems)=="r"))]
                
                reactivevalue$scores <- pems
                
                updateSelectizeInput(inputId = "gene_name",
                                     choices = rownames(reactivevalue$se),
                                     selected = NULL)
                
                pems <- as.data.table(pems)
                
                rownames(pems) <- rownames(reactivevalue$se)
                
                output$pem <- renderDT(pems
                )
                
                output$download_pems <- downloadHandler(
                    filename = function() {
                        paste("PEM_scores", ".csv", sep = "")
                    },
                    content = function(file) {
                        write.csv(reactivevalue$scores,file)
                    }
                )
            })
            
        }
        ) 
        
        observeEvent(input$plot, {
            
            cols <- c(brewer.pal(length(colnames(reactivevalue$scores)),name = 'Pastel1'))
            
            if (input$plot_type == 'scatter') {
                
                x <- round(length(colnames(reactivevalue$scores))/2)
                y <- round(length(colnames(reactivevalue$scores)))
                
                output$plots <- renderPlot({
                    
                    par(mfrow=c(x,y)) 
                    
                    for (i in 1:length(colnames(reactivevalue$scores))) {
                        
                        counts <- log10(reactivevalue$data[,i])
                        
                        pems <- reactivevalue$scores[,i] 
                        smoothScatter(counts,
                                      pems,
                                      ylab="PEM Scores", 
                                      xlab="log10 Counts",
                                      main=colnames(reactivevalue$data[i]))
                        points(log10(reactivevalue$data[input$gene_name,i]),
                               reactivevalue$scores[input$gene_name,i],
                               col="red",
                               pch=16,
                               cex=0.5)
                    }
                }
                )
                
            }
            
            if (input$plot_type == 'bar') {
                
                output$plots <- renderPlot(barplot(reactivevalue$scores[input$gene_name,],
                                                   xlab='Sample Condition',
                                                   ylab='PEM Score',
                                                   main=input$gene_name,
                                                   col=cols)
                )
            }
            
            zero <- c()
            negative <- c()
            positive <- c()
            for (i in colnames(reactivevalue$scores)) {
                
                total <- length(reactivevalue$scores[,i])
                
                zero <- c(zero,round(length(which(reactivevalue$scores[,i]==0))/total,digits = 2))
                negative <- c(negative,round(length(which(reactivevalue$scores[,i]<0))/total,digits = 2))
                positive <- c(positive,round(length(which(reactivevalue$scores[,i]>0))/total,digits = 2))
            }
            
            pem_percents <- abind(zero,negative,positive,along=2)
            
            rownames(pem_percents) <- colnames(reactivevalue$scores)
            
            colnames(pem_percents) <- c('% zero','% negative','% positive')
            
            output$tab <- renderDT(pem_percents)
            
        }
        )
        
        
    })
