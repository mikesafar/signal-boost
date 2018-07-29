#######################################################################
# Mike Safar - CorpusSummary
# ----------------------------------------------
# Copyright (C) 2018.  All Rights Reserved.  
########################################################################

library(R6)
library(logging)
library(tm)
library(scales)
library(dplyr)

CorpusSummary <- R6Class("CorpusSummary",
      
      public <- list(
        
    #constructor
        initialize = function(corpus, 
                              pre.stemmed.corpus = NULL, 
                              weight.terms = T,
                              sparse.maximal = 0.5,
                              method = "euclidian",
                              k.clusters = 5, 
                              k.rounds = 20, 
                              min.words.per.doc = NULL,
                              ...) {
            private$corpus <- corpus
            private$weight_DTM <- weight.terms
            private$method <- method
            private$pre_stemmed_corpus <- pre.stemmed.corpus
            private$sparse_maximal <- sparse.maximal
            private$kClusters <- k.clusters
            private$kRounds <- k.rounds
            private$min_words_per_doc <- min.words.per.doc
      
            # Not sure how to handle this, but seems like i have to ensure it's there
            addHandler(writeToConsole)
                  
            private$processCorpus(corpus, ...)  
        },
    
    #getCorpus
        getCorpus = function() {private$corpus},
    
    #getDtm
        getDTM = function(sparse = TRUE) {
            if (sparse)  
                private$dtm_sparse
            else
                private$dtm
        },
    
    #getMostFrequentTerms
        getMostFrequentTerms = function(sparse = TRUE) {
            private$getMostFreqTerms(self$getDTM(sparse))
        },
    
    #getKmeansResults
        getKmeansResults = function() {private$kResults},
    
    #getDist
        getDist = function(completions = FALSE) {
            dist <- private$dist
            
            if (completions & !is.null(private$stem_completions)) {
                ct <- private$stem_completions
                labels <- attr(dist, "Labels")
                new.labels <- as.vector(ct[ct$stem %in% labels,]$completion)
                attr(dist, "Labels") <- new.labels
            }
            
            return(dist)
        },
        
    #getSummary    
        getSummary = function(min.words.per.doc = NULL) {
            private$summary
        },
    
        getStemCompletionTable = function() {private$stem_completions},
    
        getTermsFromDoc = function(doc.number) {
            if (is.na(doc.number) || !is.numeric(doc.number) || is.null(doc.number) || identical(doc.number, integer(0)) || !(doc.number > 0))
                return(NULL)
          
            dtm <- self$getDTM(sparse = FALSE)
            
            #assumes tfidf
            pruned.dtm <- dtm[doc.number,as.vector(dtm[doc.number,] > 0)]
            weighted.terms <- private$getMostFreqTerms(pruned.dtm)
            
            if (!is.null(private$stem_completions)) {
              weighted.terms <- merge(weighted.terms, private$stem_completions, by.x = "word", by.y = "stem")
              rownames(weighted.terms) <- weighted.terms$word
              weighted.terms$word <- weighted.terms$completion
              weighted.terms$completion <- NULL
            } else {
              rownames(weighted.terms) <- weighted.terms$word
            }

            weighted.terms <- weighted.terms[order(-weighted.terms$freq),]
                        
            return(weighted.terms)
        }
    
      ),
  
      private <- list(
        
    #summary    
        corpus = "SimpleCorpus",
        pre_stemmed_corpus = "SimpleCorpus",
        dtm = "DocumentTermMatrix",
        dtm_sparse = "DocumentTermMatrix",
        dist = "dist",
        method = "character",
        weight_DTM = "logical",
        word_frequencies = NULL,
        stem_completions = NULL,
        sparse_maximal = "numeric",
        kClusters = "numeric",
        kRounds = "numeric",
        kResults = "kmeans",
        min_words_per_doc = "numeric",
        summary = "list",
        summaryCompleted = FALSE,
        
    #processCorpus 
        processCorpus = function(corpus, ...) {
            
            loginfo("PROCESSING CORPUS WTIH DOCUMENTS: %d", length(corpus))
          
            loginfo("...creating Matrix...")
            private$dtm <- DocumentTermMatrix(corpus, ...)
            
            loginfo("......found %s terms...", comma_format()(length(private$dtm$dimnames$Terms)))
            
            loginfo("...getting stem completions...") 
            private$stem_completions <- private$getCompletionTable()
            
            loginfo("...weighting the matrix...")
            if (private$weight_DTM) {
                suppressWarnings(private$dtm <- weightTfIdf(private$dtm))  
            }
            
            loginfo("...removing sparse terms at maximal of %f...", private$sparse_maximal)
            private$dtm_sparse <- removeSparseTerms(private$dtm, private$sparse_maximal)
           
            loginfo("......reduced to %s terms...", comma_format()(length(private$dtm_sparse$dimnames$Terms)))
            
            loginfo("...creating %s distance matrix...", private$method)
            private$dist <- dist(t(private$dtm_sparse), private$method)
            
            loginfo("...finding %d kmeans clusters over %d rounds...", private$kClusters, private$kRounds)
            private$kResults <- kmeans(private$dist, private$kClusters, private$kRounds)
            
            loginfo("...composing summary...")
            private$composeSummary(private$min_words_per_doc)
            
            loginfo("...DONE processing.")
        },
    
    #getMostFreqTerms
        getMostFreqTerms = function(dtm = NULL, wordList = NULL) {
            if (is.null(dtm)) 
                m <- private$dtm_sparse
            else
                m <- dtm

            ft.words <- colnames(m)
            ft <- data.frame(word=colnames(m), freq=col_sums(m), row.names=ft.words, stringsAsFactors = F)

            if (!is.null(wordList) && length(wordList) > 0) {
                ft <- ft[which(ft$word %in% wordList),]
            }
              
            if (!is.null(private$stem_completions)) {
                ft <- private$getCompletedWords(ft)
            }
            
            ft <- ft[order(ft$freq, decreasing=TRUE),]  

            return(ft)
        },
    
    #getCompletedWords -- WILL RETURN FREQUENCIES FROM PRE-STEMMED CORPUS
        getCompletedWords = function(stems) {
            completions = private$stem_completions
            stopifnot(!is.null(completions))
            
            cw <- merge(x=stems, y=completions, by.x="word", by.y="stem")
            rownames(cw) <- cw$word
            cw$word <- cw$completion
            cw$completion <- NULL
            
            return(cw)
        },
      
    #getCompletionTable
        getCompletionTable = function(corpus = NULL) {
            c <- NULL
            if (is.null(corpus))
                c <- private$pre_stemmed_corpus
            else
                c <- corpus
            
            if (!is.null(c)) {
                originals <- DocumentTermMatrix(c)
                completion <- col_sums(originals) %>% sort(decreasing=T)
                stem <- stemDocument(names(completion), language = "en")
                completionTable <- as.data.frame(cbind(stem, names(completion)), row.names=stem, stringsAsFactors = T)
                colnames(completionTable) <- c("stem", "completion")
                completionTable <- completionTable[!duplicated(completionTable$stem),]
                
                return(completionTable)
            }
          
            return(NULL)
        },
    
    #composeSummary    
        composeSummary = function(min.words.per.doc = NULL) {
            stopifnot(!is.null(min.words.per.doc) || !is.numeric(min.words.per.doc))
          
            results <- private$kResults
            
            clusterSummaries <- list()
            for (i in 1:private$kClusters) {
                clusterSummaries[[i]] <- private$composeClusterSummary(i, min.words.per.doc)
            }
            
            private$summary <- clusterSummaries
            private$summaryCompleted <- TRUE
            
            return(private$summary)
        },
        
    #composeClusterSummary
        composeClusterSummary = function(clusterNumber, min.words.per.doc = NULL) {
            stopifnot(!is.null(min.words.per.doc) || !is.numeric(min.words.per.doc))
          
            results <- private$kResults
            dtm <- as.matrix(private$dtm_sparse)
            rownames(dtm) <- 1:nrow(dtm)
          #order of operations important here!
            termList <- names(which(results$cluster == clusterNumber))
            termFreqTable <- as.data.frame(private$getMostFreqTerms(dtm, termList))
            
            min <- min.words.per.doc
            if (is.null(min)) min <- 0
                
            relevant.docs <- which(rowSums(as.matrix(dtm[,termList]) > 0) > min)
          
            docList <- rowSums(as.matrix(dtm[relevant.docs,termList])>0)
            docList <- as.integer(names(which(docList[order(docList, decreasing=TRUE)] > 0)))
            
            loginfo("--- *** Cluster %d: %s Documents, %d Terms", clusterNumber, comma_format()(length(docList)), length(termFreqTable$word))
                    
            return(
                list(
                    "docList" = docList,
                    "termList" = termFreqTable
                )
            )
        }
    
      )
  )
    


