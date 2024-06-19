#' Merge dataframes by shared rows
#' Dataframes to be merged need to have rownames with the TFs or other entities being merged on.
#' @param df1 Dataframe 1, essentially a named vector in dataframe format
#' @param df2 Dataframe 2, essentially a named vector in dataframe format with
#' @param margin The direction of merging dataframes only implemented for rows
#' @export
vector_smartmerge <- function(df1,df2,margin="row"){
    if(margin == 'row'){
        row_names <- union(rownames(df1),rownames(df2))
        df_int <- cbind(df1,0)
        col_index <- dim(df_int)[2]
        df_int[rownames(df2),col_index] <- df2[,1]
        colnames(df_int)[col_index] <- colnames(df2)
    }
    else{print("Haven't worked this out yet")}
             return(df_int)}


#' Calculate Metric
#' Calculates 3 major network features from edgelists of TFs and gene targets in igraph.rds format or txt files. Edgelists need to be in tab separated format with TF in column 1, Gene in column 2. Additional columns are ignored in the current implementation.
#' @param indir Directory containing edgelists for all networks. MUST BE IN .rds for igraph objects or tab separated .txt files. It helps if the the filenames include the sample ID followed by the file extension.
#' @param outfile File to write final metric matrix
#' @param metric Metric to calculate
#' @param in_ext File extension of edgelists
#' @param out_ext 
#'
#' @export 

calc_metric <- function(indir, outfile,metric=c("outdegree","indegree","betweenness"),in_ext="_graph.rds",out_ext=".rds"){
##    require(igraph)

    file_list <- list.files(indir,pattern=in_ext)

    metric_df <- data.frame(0)
    rownames(metric_df) <- "Empty"
    str(file_list)

    print(paste0("Calculating ",metric))
    
### Checking file formats for compatibility
    
    check <- ifelse(length(grep("\\.rds",file_list))>0,"rds",ifelse(length(grep("\\.txt",file_list))>0,"txt","unknown"))



    for(i in file_list){

        if(check=="rds") {graph <- readRDS(paste0(indir,i))} else if (check=="txt"){ graph <- read.table(paste0(indir,i),header=F, sep='\t')
                                                                                     graph <- graph_from_edgelist(as.matrix(graph[,c(1,2)]))} else if (check=="unknown") { print("This is an unsupported file extension format")}                                                                                   
        if(metric=="outdegree"){
            metric_out <- degree(graph,mode='out')
        } else if(metric=="indegree"){
            metric_out <- degree(graph,mode="in")
        }
            else if ( metric=="betweenness"){ metric_out <- betweenness(graph)
                                          }
    metric_out <- metric_out[which(metric_out >0)]
        metric_out <- as.data.frame(metric_out, row.names=names(metric_out))
        colnames(metric_out) <- i

  

    metric_df <- vector_smartmerge(metric_df,metric_out,'row')

    print(paste0("Finished processing ", grep(i,file_list)," of ", length(file_list)))
    }
    colnames(metric_df) <- gsub(paste0("_Final|_PIQ|",ext),"",colnames(metric_df))
    metric_df <- metric_df[-1,-1]
    metric_df[is.na(metric_df)] <- 0

saveRDS(metric_df,outfile)
}

get_regulator_TFs <- function(graph,Gene){
   ## require(igraph)
    index <- which(names(V(graph)) == Gene)
    if(length(index) >0 ){    regulators <- names(which(graph[,index]!=0))} else{ regulators <- ""}
    return(regulators)}


#' Get Target Genes
#' Identifies target genes of a given TF
#' @param graph igraph object for a given sample/cell line
#' @param TF TF of interest
#'
#' @export 

get_regulatory_targets <- function(graph,TF){
   ## require(igraph)
    index <- which(names(V(graph)) == TF)
    if(length(index) > 0){    targets <- names(which(graph[index,]!=0))} else {targets <- ""}
    return(targets)}


#' Calculate Metric Stability
#' Calculates 3 major network features from edgelists of TFs and gene targets in igraph.rds format or txt files. Edgelists need to be in tab separated format with TF in column 1, Gene in column 2. Additional columns are ignored in the current implementation.
#' @param indir Directory containing edgelists for all networks. MUST BE IN .rds for igraph objects or tab separated .txt files. It helps if the the filenames include the sample ID followed by the file extension.
#' @param outfile File to write final metric matrix. Can be NULL to return matrix in session
#' @param fraction Fraction of edges in a given network to remove for evaluating stability
#' @param metric Metric to calculate
#' @param verbose Suppress warnings and output logging
#' @param in_ext File extension of edgelists
#' @param out_ext File extension of output file
#'
#' @export 


calc_metric_stability <- function(indir, outfile=NULL,fraction,metric=c("outdegree","indegree","betweenness"),verbose=FALSE,in_ext="_graph.rds",out_ext=".rds"){
    require(igraph)

    file_list <- list.files(indir,in_ext)

    metric_df <- data.frame(0)
    rownames(metric_df) <- "Empty"

    print(paste0("Calculating ",metric))
    
### Checking file formats for compatibility
    
    check <- ifelse(length(grep("\\.rds",file_list))>0,"rds",ifelse(length(grep("\\.txt",file_list))>0,"txt","unknown"))



    for(i in file_list){

        if(check=="rds") {graph <- readRDS(paste0(indir,i))} else if (check=="txt"){ graph <- read.table(paste0(indir,i),header=F, sep='\t')
                                                                                     graph <- graph_from_edgelist(as.matrix(graph[,c(1,2)]))} else if (check=="unknown") { print("This is an unsupported file extension format")}
        
        graph <- delete_edges(graph,sample(V(graph),(length(V(graph))*fraction)))

        if(metric=="outdegree"){
            metric_out <- degree(graph,mode='out')
        } else if(metric=="indegree"){
            metric_out <- degree(graph,mode="in")
        }
            else if ( metric=="betweenness"){ metric_out <- betweenness(graph)
                                          }
    metric_out <- metric_out[which(metric_out >0)]
        metric_out <- as.data.frame(metric_out, row.names=names(metric_out))
        colnames(metric_out) <- i

    #str(metric_out)

    metric_df <- vector_smartmerge(metric_df,metric_out,'row')
    if(verbose== TRUE){
    print(paste0("Finished processing ", grep(i,file_list)," of ", length(file_list)))}
    }
    colnames(metric_df) <- gsub(paste0("_Final|_PIQ|",in_ext),"",colnames(metric_df))
    metric_df <- metric_df[-1,-1]

if(is.null(outfile)){ return(metric_df)} else { if (out_ext ==".rds") {saveRDS(metric_df,outfile)} else { write.table(metric_df, outfile, quote=F, sep='\t')}}
}

#' Calculate TF-Score
#' Calculates TF Score for a set of matrices where TFs are rows and samples/cell lines are columns. This works best for RANKED matrices. Only keeps TFs present in all provided matrices.
#' @param outdegree_df Matrix of outdegree centrality for TFs
#' @param expression_df Matrix of expression for TFs
#' @param betweenness_df Matrix of betweenness centrality for TFs is optional
#' @param method Simple average of the matrices or average + standard deviation. 
#' @export 

calc_TF_score <- function(outdegree_df,expression_df, betweenness_df=NULL,method=c("simple","complex")){
    if(is.null(betweenness_df)){

    common_genes <- intersect(rownames(outdegree_df),rownames(expression_df))
    colname_index <- intersect(colnames(outdegree_df),colnames(expression_df))

    outdegree_df <- outdegree_df[common_genes,colname_index]
    expression_df <- expression_df[common_genes,colname_index]

    if(method == "simple"){
    out_df <- (outdegree_df+expression_df)/2
    out_df[is.na(out_df)] <- 0} else if(method=="complex"){
        out_df <- (outdegree_df+expression_df)/2
        int_df <- matrix(mapply(function(x,y) sd(c(x,y)),unlist(outdegree_df),unlist(expression_df)),ncol=ncol(outdegree_df),nrow=nrow(outdegree_df))

        out_df <- out_df+int_df

        #out_df <- rank_mat(out_df)
        out_df[is.na(out_df)] <- 0
        out_df <- as.data.frame(out_df)}

}
    else{
        common_genes <- Reduce(intersect, list(rownames(outdegree_df), rownames(expression_df), rownames(betweenness_df)))
        colname_index <- Reduce(intersect, list(colnames(outdegree_df), colnames(expression_df), colnames(betweenness_df)))
        outdegree_df <- outdegree_df[common_genes,colname_index]
        expression_df <- expression_df[common_genes,colname_index]
        betweenness_df <- betweenness_df[common_genes,colname_index]

    if(method == "simple"){
    out_df <- (outdegree_df+expression_df+betweenness_df)/3
    out_df[is.na(out_df)] <- 0} else if(method=="complex"){
        out_df <- (outdegree_df+expression_df+betweenness_df)/3
        int_df <- matrix(mapply(function(x,y,z) sd(c(x,y,z)),unlist(outdegree_df),unlist(expression_df),unlist(betweenness_df)),ncol=ncol(outdegree_df),nrow=nrow(outdegree_df))

        out_df <- out_df+int_df

        #out_df <- rank_mat(out_df)
        out_df[is.na(out_df)] <- 0
        out_df <- as.data.frame(out_df)}





    }

    return(out_df)}

#' Rank Matrices
#' Ranks a matrix by either row or column such that high values in original matrix are given rank 1-n
#' @param matrix Matrix of whatever you want to rank, must be numeric
#' @param margin How to rank matrix, row calculates the rank within rows, column calculates the rank within columns
#' @export 

rank_mat <- function(matrix,margin="column"){
    if(any(is.na(matrix)) == TRUE){
        matrix[is.na(matrix)] <- 0
        print("There are NA values in your matrix, converting to 0")} else{ matrix= matrix}

    if(margin == "column"){
        mat <- apply(matrix,2, function(x) (length(x)+1)-rank(x))}
    else {         mat <- t(apply(matrix,1, function(x) (length(x)+1)-rank(x)))}
    mat <- as.data.frame(mat)
    mat[is.na(mat)] <- 0
    return(mat)}
