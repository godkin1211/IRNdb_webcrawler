#!/usr/bin/env Rscript --vanilla
###################################################################
# Name: web_crawler.r                                             #
# Description: This script is used to scrape lncRNA-mRNA target-  #
#              informations from IRNdb.                           #
# Required packages: rvest, dplyr, RCurl, jsonlite, V8            #
# Author: Chiachun Chiu <chiachun.chiu@gmail.com>                 #
#                       <jiajyun.ciou@gmail.com>                  #
# Date: 2017-11-8                                                 #
###################################################################

# Importing required pkgs
library(rvest)
library(RCurl)
library(rvest)
library(jsonlite)
library(V8)

# Settings
mainDir <- file.path("~/Projects/Henan")
lncRNA_table_file <- file.path("~/Desktop/irndb_export.csv") # Necessary file
target_gene_table_header <- c("Symbol", "Name", "Species", "Source", "TargetSource", "PMID")

# Self-defined function
partitionVec <- function(vec, partitionSize=10) {
    vec.len <- length(vec)
    if (vec.len <= 20) stop("Empty or small vector!")
    vec.parted <- split(vec, as.integer((seq(vec.len) - 1)/partitionSize))
    return(vec.parted)
}

web_crawler <- function(links) {
    if (!is.vector(links)) stop("Type Error! Please check if given links are correct!")
    if (length(links)==0 | links=="") stop("Empty vector! Please check if given links are correct!")
    if (!grepl("^http", links[1])) stop("The input has to be a valid link beginning with 'http'!")
    htmltxt.vec <- c()
    # If the number of links is larger than 20, divide and conquer!
    if (length(links) <= 20) {
        htmltxt.vec <- getURIAsynchronous(links)
    } else {
        links.parted <- partitionVec(links, 10)
        for (i in seq(length(links.parted))) {
            new.htmltxt <- getURIAsynchronous(links.parted[[i]])
            htmltxt.vec <- c(htmltxt.vec, new.htmltxt)
            Sys.sleep(10)
        }
    }
    if (length(htmltxt.vec) != length(links)) stop("Some pages missed!")
    cat("Web-scraping job has been completed!\n")
    return(htmltxt.vec)
}

html2table <- function(htmlTextObj, lncRNA, target_gene_table_header) {
    if (is.null(lncRNA) | lncRNA=="") stop("Please check the variable is defined correctly: lncRNA.")
    if (is.null(target_gene_table_header) | target_gene_table_header=="") stop("Please check the variable is defined correctly: target_gene_table_header")
    if (is.null(htmlTextObj) | htmlTextObj=="") stop("Please check the variable is defined correctly: htmlTextObj")
    htmlcontent <- read_html(htmlTextObj) %>% html_nodes("div.dataTable_wrapper script")
    ct <- v8()
    ct$eval(gsub('\\\n   ', '', gsub('\\\n    \\$.*','', htmlObj)))
    target_gene_table <- as.data.frame(ct$get("dataSetT"))
    #target_gene_table <- htmlcontent %>%  
    #                     html_nodes("script:contains('var dataSetT')") %>% 
    #                     gsub(".*var dataSetT = ", "", .) %>% 
    #                     gsub("a>\\']];.*\\$\\(document\\).*", "a>\\']]", .) %>% 
    #                     gsub("\'", "\"", .) %>% 
    #                     gsub("=\"","=\'",.) %>% 
    #                     gsub("\" ", "\' ", .) %>% 
    #                     gsub("\">", "\'>", .) %>%
    #                     gsub("5\\\\\' ", "", .) %>% 
    #                     gsub("3\\\\\'", "", .) %>% 
    #                     fromJSON %>% 
    #                     as.data.frame
    colnames(target_gene_table) <- target_gene_table_header
    target_gene_table$lncRNA <- lncRNA
    return(target_gene_table)
}

# Main
cat("###############################################################\n")
cat("# Script start!                                               #\n")
cat("###############################################################\n\n")
cat("* Reading lncRNA_list.......\n")
lncRNA_table <- read.csv(lncRNA_table_file, stringsAsFactors = FALSE)
lncRNA_list <- lncRNA_table$Symbol
if (!is.vector(lncRNA_list)) lncRNA_list <- as.character(lncRNA_list)
cat("* Generating links for web-scraping....\n")
irndb_links <- paste0("http://compbio.massey.ac.nz/apps/irndb/lncrna/",lncRNA_list,"?type=t")
cat("* Starting to crawl target gene information from IRNdb......")
html_txt_content_vec <- web_crawler(irndb_links)
target_gene_table <- NULL
pb <- txtProgressBar(min=1, max=length(html_txt_content_vec), style=3)
cat("* Starting to parse the html-contents we scrapped.....\n")
for (i in 1:length(html_txt_content_vec)) {
    lncRNA <- lncRNA_list[i]
    #cat("  --- processing the", i, "lncRNA -", lncRNA, ".....\n")
    target_page <- html_txt_content_vec[i]
    output_table <- html2table(target_page, lncRNA, target_gene_table_header)
    target_gene_table <- rbind(target_gene_table, output_table)
    setTxtProgressBar(pb, i)
}
close(pb)
cat("Page-parsing job has been completed!\n")
cat("* Writing output file........\n")
output_file <- "IRNdb_lncRNA_target_gene_table.txt"
write.table(target_gene_table, file.path(mainDir, output_file), row.names = FALSE, sep="\t")
cat("Total jobs finished! Goodbye!\n")