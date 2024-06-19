ANF_qqplot <- function (pval_df,TF_list, regulating_TFs,Gene){
    pval_df$observed <- -log10(pval_df[,1])
    pval_df$expected <- -log10(runif(lengt(pval_df[,1])))
    pval_df$TF <- "No"
    pval_df$ID <- rownames(pval_df)
    TFs <- TF_list
    pval_df$TF[which(rownames(pval_df) %in% TFs)] <- "Yes"
    pval_df$TF[which(rownames(pval_df) %in% regulating_TFs)] <- paste0(Gene,"\nregulating")
    qplot <- ggplot(pval_df, aes(sample=observed))+stat_qq()
    qplot_df <- ggplot_build(qplot)$data[[1]]
    qplot_df$TF <- pval_df$TF[order(pval_df$observed)]
    qplot_df$ID <- pval_df$ID[order(pval_df$observed)]
    qplot_df$theoretical <- sort(-log10(runif(length(qplot_df$theoretical))))
    p <- ggplot(qplot_df,aes(theoretical,sample,label=ID))+geom_point(aes(color=qplot_df$TF))+geom_text(nudge_x=-0.2,size=1.5,aes(color=qplot_df$TF),check_overlap=TRUE)+geom_abline(linetype="dotted",slope=1,intercept=0)+geom_hline(aes(yintercept=-log10(0.05/length(qplot_df$TF))), color="coral", size=2)

    return(p)}

ascat_to_granges <- function(ascat_summary){ ascat_summary <- read.table(ascat_summary,header=T,stringsAsFactors=F,sep=',');ascat_summary[,1] <- gsub(23,"X",ascat_summary[,1]);ascat_summary[,1] <- gsub(24,"Y",ascat_summary[,1]);gr <- GRanges(seqnames=Rle(ascat_summary[,1]), ranges= IRanges(start=ascat_summary[,2], end=ascat_summary[,3]),strand="*",CN=ascat_summary[,6], Sample_ID=ascat_summary[,8]); gr <- gr.chr(gr);                               return(gr)}

bed_to_granges <- function(bed,header=FALSE){
    require(GenomicRanges)
    bed <- read.delim(bed,header=header, stringsAsFactors=F, sep="\t")

                                 gr <- GRanges(seqnames=Rle(bed[,1]), ranges= IRanges(start=bed[,2], end=bed[,3]),strand="*")
                                 return(gr)}

bed_to_granges_enhancer <- function(bed){ require(GenomicRanges)
                                          bed <- read.table(bed, stringsAsFactors=F, sep="\t")
                                          gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",Target=bed$V4) ;return(gr)}
bed_to_granges_enhancer_extra <- function(bed){ require(GenomicRanges)
                                                bed <- read.table(bed, stringsAsFactors=F, sep="\t"); gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",Target=bed$V4, sample=bed$V5) ;return(gr)}


bed_to_granges_extra <- function(bed){ require(GenomicRanges)
                                       bed <- read.table(bed, stringsAsFactors=F, sep="\t")
                                 gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",sample=bed$V7) ;return(gr)}
bed_to_granges_SNV_extra <- function(bed){ require(GenomicRanges)
                                                bed <- read.table(bed, stringsAsFactors=F, sep="\t"); gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",Ref=bed$V4, Alt=bed$V5, Sample=bed$V6) ;return(gr)}


chrom_hmm_to_granges <- function(bed) {
    require(GenomicRanges)
    bed <- read.table(bed, stringsAsFactors=F, sep="\t")
    bed$V5 <- as.numeric(unlist(lapply(bed$V4, function(x) unlist(strsplit(x,"_"))[1])))
    gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),strand="*",State=bed$V4,score=bed$V5) ;return(gr)}

gr_to_bed <- function(gr,outfile=NULL, metadata=FALSE,verbose=FALSE,header=FALSE,biostrings=FALSE,append=FALSE,drop_strand=TRUE){
    require(GenomicRanges)

if(length(gr) >0){
    Chrom <- as.character(gr@seqnames)
    Start <- gr@ranges@start
    End <- gr@ranges@start+gr@ranges@width
    Strand <- gr@strand

                                   if(metadata==FALSE){
                                       df <- data.frame(Chrom,Start,End)}
                                   else if(metadata==TRUE){
                                       metadata <- as.data.frame(gr@elementMetadata@listData,stringsAsFactors=F)
                                       if(biostrings==TRUE){
                                           metadata <- metadata[,setdiff(colnames(metadata),c("ALT.group","ALT.group_name"))]}
                                       df <- data.frame(Chrom,Start,End,Strand,metadata,stringsAsFactors=F)

                                       }
    if(is.null(outfile)){
        return(df)} else{
            if(drop_strand== TRUE){ df["Strand"] <- NULL} else{ df <- df}
            write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=header,append=append)}
                                   if(verbose== TRUE) {print(head(df))} else{}} else{ print("Zero length gr. Stopping")}}

links_to_bed <- function(link_names, peak_file,outfile, metadata=FALSE){
    require(GenomicRanges)
    split1 <- link_names
    split1 <- unlist(strsplit(split1,"~"))
    peaks <- split1[seq(1,length(split1),2)]
    targets <- split1[seq(2,length(split1),2)]
    index <- which(peaks %in% names(peak_file))
    true_peaks <- peaks[index]
    true_targets <- targets[index]
    #str(true_peaks)
    out_gr <- peak_file[true_peaks]
    Chrom <- as.character(out_gr@seqnames)
    Start <- out_gr@ranges@start
    End <- out_gr@ranges@start+out_gr@ranges@width
    if(metadata==FALSE){
        df <- data.frame(Chrom,Start,End)}
    else if(metadata==TRUE){
        metadata <- as.data.frame(out_gr@elementMetadata@listData)
        df <- data.frame(Chrom,Start,End,metadata)
    }
    df$Identifier <- true_peaks
    df$Target <- true_targets
    df <- df[order(df[,1],df[,2]),]
    write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=F)
#    print(head(df))
}



mutations_gr_to_bed <- function(gr,outfile){
    require(GenomicRanges)
    if(length(gr) >0){
    Chrom <- as.character(gr@seqnames)
    Start <- gr@ranges@start
    End <- gr@ranges@start+gr@ranges@width
    Ref <- gr$Ref
    Alt <- gr$Alt
    Gene <- gr$Gene
    df <- data.frame(Chrom,Start,End,Ref,Alt)

    write.table(df, outfile,quote=F, sep='\t', row.names=F,col.names=F)} else{ print("Zero length gr")}}



Binarize_CN_table_subtype <- function(CN_table){
    CN_table$Binarized_Subtype <- 0
    CN_table[CN_table[,1] >6,"Binarized_Subtype"] <- 1
    return(CN_table)}


CN_to_granges <- function(CN_df){ require(GenomicRanges)
                                  CN_df <- read.table(CN_df,header=F,stringsAsFactors=F,sep='\t');CN_df[,1] <- gsub(23,"X",CN_df[,1]);CN_df[,1] <- gsub(24,"Y",CN_df[,1]); CN_df[,6] <- gsub("a|a_2|b|c","",CN_df[,6]);gr <- GRanges(seqnames=Rle(CN_df[,1]), ranges= IRanges(start=CN_df[,2], end=CN_df[,3]),strand="*",CN=CN_df[,5], Sample_ID=CN_df[,6]); gr <- gr.chr(gr);                               return(gr)}



Copy_Number_Eval <- function(CN_gr,interval_gr,expression_matrix,Gene){
    CNA_table <- as.data.frame(expression_matrix[,Gene])
    samples_intersecting_Gene <- intersect_with_metadata(CN_gr,interval_gr)
    seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)

    rownames(CNA_table) <- rownames(expression_matrix)
    CNA_table$Copy_Number <- NA
    colnames(CNA_table) <- c(Gene,"Copy_Number")
    CNA_table$ID <-  rownames(CNA_table)
    for(i in rownames(CNA_table)){ index <- which(samples_intersecting_Gene$Sample_ID == i); CN <- mean(samples_intersecting_Gene$CN[index]);CNA_table[i,]$Copy_Number <- CN}
    CNA_table$Subtype <- NA; for(i in rownames(CNA_table)){sample_subtype <- expression_matrix[i,"Status"];CNA_table[i,]$Subtype <- sample_subtype}
return(CNA_table)}

Cross_validation <- function(matrix,num_of_folds,featurelist,Gene){
    fold_size <- round(dim(matrix)[1]/num_of_folds)
    out_list <- list()
    blacklist <- c()

    for(i in 1:num_of_folds){
        test <- sample(setdiff(rownames(matrix),blacklist),fold_size)
        blacklist <- c(blacklist,test)
        train <- setdiff(rownames(matrix),test)
        assign(paste0("model_",i),multiple_regression(matrix[train,],featurelist,Gene))
        test_result <- predict.glm(get(paste0("model_",i)),matrix[test,])
        df_out <- data.frame(matrix[test,Gene],test_result)
        colnames(df_out) <- c("Actual","Prediction")
        df_out$Features <- paste0(featurelist,collapse=",")
        out_list <- append(out_list,list(df_out))
        print(paste0("Finished Processing Fold ",i))
    }
    return(out_list)}


Condense_cross_validation_output <- function(Cross_validation_output,feature_info){
    final_df <- data.frame()
    cor_val_vec <- c()
    for(i in 1:length(Cross_validation_output)){
        sub <- Cross_validation_output[[i]]
        cor_val <- cor(sub[,1],sub[,2],use='complete.obs')^2
        Fold <- paste0("Fold_",i)

        out_df <- data.frame(cor_val,Fold,unique(sub$Features))
        final_df <- rbind(final_df,out_df)
    }
    colnames(final_df) <- c("R_Squared","Fold","Features")
    final_df$Feature_Type <- feature_info
   return(final_df)
}



Plot_cross_validation_output <- function(Cross_validation_output,featurelist){
    final_df <- data.frame()
    cor_val_vec <- c()
    for(i in 1:length(Cross_validation_output)){
        sub <- Cross_validation_output[[i]]
        cor_val <- cor(sub[,1],sub[,2],use='complete.obs')^2
        Fold <- paste0("Fold_",i)

        out_df <- data.frame(cor_val,Fold,sub$Features)
        final_df <- rbind(final_df,out_df)
    }
    colnames(final_df) <- c("R_Squared","Fold","Features")
    final_df$Features <- paste0(featurelist,collapse="-")
    p <- ggplot(final_df,aes_string(x="Fold",y="R_Squared"))+geom_point()+geom_hline(yintercept=mean(final_df[,1]),size=2, alpha=0.5,color="coral", linetype=2)+labs(title=paste0("Model Performance in ",length(Cross_validation_output),"-fold cross-validation"),subtitle=final_df$Features)
    print(p)
    dev.off()

}

Plot_condensed_cross_validation_output <- function(condensed_output,fig_title){
    f_vec <- c()
    for(i in 1:length(condensed_output$Features)){
        f_len <- length(unlist(strsplit(as.character(condensed_output$Features[i]),",")))
        f_vec <- c(f_vec,f_len)}
    condensed_output$Model <- as.factor(f_vec)

    p <- ggplot(condensed_output, aes(x=Model,y=R_Squared,fill=Feature_Type))+geom_violin(trim=TRUE,draw_quantiles=c(0.25,0.5,0.75))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(title=fig_title,x="Number of model features")+geom_smooth(data=condensed_output,aes(x=as.numeric(Model),y=R_Squared))
    print(p)
    dev.off()
}

get_targeting_enhancers <- function(gr,metadata_column,gene){ require(GenomicRanges)
                                                              out <- gr[which(gr@elementMetadata[,metadata_column] == gene)]; return(out)}

ID_enhancer_mutations <- function(enhancer_gr, Mutation_calls){
    require(GenomicRanges)
    infile <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    infile_gr <- GRanges(seqnames=Rle(c(infile$Chrom)), ranges=IRanges(start=infile$start, end=infile$end), Gene=infile$Gene, Ref=infile$Ref, Alt=infile$Alt, Sample=infile$Sample)
    mutations_gr <- infile_gr[findOverlaps(infile_gr,enhancer_gr)@from]
    return(mutations_gr)
}

ID_gene_mutations <- function(gene_of_interest, Mutation_calls){
    require(GenomicRanges)
    infile <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    infile_gr <- GRanges(seqnames=Rle(c(infile$Chrom)), ranges=IRanges(start=infile$start, end=infile$end), Gene=infile$Gene, Ref=infile$Ref, Alt=infile$Alt, Sample=infile$Sample)
    mutations_gr <- infile_gr[which(infile_gr$Gene == gene_of_interest),]
    return(mutations_gr)
}

intersect_with_metadata <- function(gr1,gr2, ignore.strand=FALSE){
    require(GenomicRanges)
    if(ignore.strand == FALSE){ out <- gr1[findOverlaps(gr1,gr2)@from,]} else if( ignore.strand ==TRUE){ out <- gr1[findOverlaps(gr1,gr2,ignore.strand=TRUE)@from,]}
    ; return(out)}

setdiff_with_metadata <- function(gr1,gr2){ require(GenomicRanges)
                                            int <- findOverlaps(gr1,gr2); out <- gr1[setdiff(1:length(gr1),unique(int@from)),]; return(out)}

Collapse_gr_keep_metadata <- function(gr){
    require(GenomicRanges)
    reduced <- reduce(gr)
    final_gr <- reduced
    for (i in 1:length(reduced)){
        overlap_enhancer = findOverlaps(reduced[i], gr)
        row_names = vector()
        for( j in 1:length(overlap_enhancer)){
            names = gr[overlap_enhancer[j]@to]$Gene
            row_names = c(names, row_names)

        }

        final_gr$nearby[i] = row_names
    }
    return(final_gr)}

make_enhancer_mutation_matrix <- function(targeting_enhancer_gr, Mutation_calls){
    require(GenomicRanges)
    infile_df <- data.frame(Mutation_calls$Chrom,Mutation_calls$start,Mutation_calls$end,Mutation_calls$Gene,Mutation_calls$Ref,Mutation_calls$Alt, Mutation_calls$Sample)
    colnames(infile_df) <- c("Chrom","start","end","Gene","Ref","Alt","Sample")
    rows <- length(unique(infile_df$Sample))
    columns <- length(targeting_enhancer_gr)
    infile_df$Chrom <- gsub("23","X",infile_df$Chrom)
    infile_df$Chrom <- gsub("24","Y",infile_df$Chrom)
    infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom)), ranges=IRanges(start=infile_df$start, end=infile_df$end), Gene=infile_df$Gene, Ref=infile_df$Ref, Alt=infile_df$Alt, Sample=infile_df$Sample)
    mutation_matrix <- matrix(data=0,nrow=rows,ncol=columns)
    rownames(mutation_matrix) <- unique(infile_df$Sample)
    colnames(mutation_matrix) <- names(targeting_enhancer_gr)
    #colnames(mutation_matrix) <- paste0(seqnames(targeting_enhancer_gr),":" ,ranges(targeting_enhancer_gr)@start,"-",ranges(targeting_enhancer_gr)@start+ranges(targeting_enhancer_gr)@width)

    seqlevelsStyle(infile_gr) <- seqlevelsStyle(targeting_enhancer_gr)


    #targeting_enhancer_gr <- reduce(targeting_enhancer_gr)
    for(i in 1:length(targeting_enhancer_gr)){
                                               mutations_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr[i])@from]
                                               if (length(mutations_gr) > 0){ out <- mutations_gr
                                                                              print(paste0(length(out)," mutations in ",i,"th enhancer"))
                                                                              mutation_matrix[out$Sample,i] <- 1}
                                               else if( length(mutations_gr) <= 0){print(paste0("No intersecting mutations in ", i,"th enhancer"))}

                                           }

return(mutation_matrix)

}

make_enhancer_sv_matrix <- function(targeting_enhancer_gr, Mutation_calls, complex_flag=FALSE,blacklist=NULL){
    require(GenomicRanges)
    infile_df <- data.frame(Mutation_calls$ID,Mutation_calls$variant_type,Mutation_calls$chr_from,Mutation_calls$chr_from_bkpt,Mutation_calls$chr_from_range,Mutation_calls$chr_from_strand,Mutation_calls$chr_to,Mutation_calls$chr_to_bkpt,Mutation_calls$chr_to_range,Mutation_calls$chr_to_strand, Mutation_calls$Score,Mutation_calls$variant_type ,stringsAsFactors=F)
    colnames(infile_df) <- c("ID","Type","Chrom1","Breakpoint1","Range1","Strand1","Chrom2","Breakpoint2","Range2","Strand2", "Score","Variant_Type")

    if(complex_flag== TRUE){
        infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom1)), ranges=IRanges(start=infile_df$Breakpoint1,end=(as.numeric(infile_df$Breakpoint1)+as.numeric(infile_df$Range1))),strand=infile_df$Strand1, Chrom2=infile_df$Chrom2, Breakpoint2=infile_df$Breakpoint2, Breakpoint_range2=infile_df$Range2,Score = infile_df$Score, strand2 = infile_df$Strand2 ,Sample=infile_df$ID, Type = infile_df$Variant_Type);    print(table(infile_gr$Variant_Type))


    }

    else {
        infile_df <- infile_df[which(infile_df$Chrom1 == infile_df$Chrom2),]
        intrachrom_var <- infile_df[which(infile_df$Variant_Type != "inversion"),]
        infile_df <- intrachrom_var
        infile_gr <- GRanges(seqnames=Rle(c(infile_df$Chrom1)), ranges=IRanges(start=infile_df$Breakpoint1, end=infile_df$Breakpoint2),strand=infile_df$Strand1, Chrom2=infile_df$Chrom2, Breakpoint2=infile_df$Breakpoint2, Breakpoint_range2=infile_df$Range2,Score = infile_df$Score, strand2 = infile_df$Strand2 ,Sample=infile_df$ID, Type = infile_df$Variant_Type)  }
    table(infile_df$Variant_Type)
    rows <- length(unique(infile_gr$Sample))
    columns <- length(targeting_enhancer_gr)
    mutation_matrix <- matrix(data=0,nrow=rows,ncol=columns)
    seqlevelsStyle(infile_gr) <- seqlevelsStyle(targeting_enhancer_gr)
    rownames(mutation_matrix) <- unique(infile_gr$Sample)
    colnames(mutation_matrix) <- names(targeting_enhancer_gr)

    if(is.null(blacklist) == FALSE){
        if(class(blacklist) == "GRanges"){
            print("Applying Blacklist")
            infile_gr <- setdiff_with_metadata(infile_gr,blacklist)
       }
                        else{ print("Blacklist needs to be in GenomicRanges format")}
    }
    else if(is.null(blacklist) == TRUE){ infile_gr <- infile_gr}
    temp_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr)@from]
    print(table(temp_gr$Type))


        for(i in 1:length(targeting_enhancer_gr)){
            mutations_gr <- infile_gr[findOverlaps(infile_gr,targeting_enhancer_gr[i])@from]


                                               if (length(mutations_gr) > 0){ out <- mutations_gr
                                                                              print(paste0(length(out)," mutations in ",i,"th enhancer"))
                                                                              pos <- out[which(out$Type != "deletion"),]
                                                                              neg <- out[which(out$Type == "deletion"),]
                                                                              mutation_matrix[neg$Sample,i] <- -1

                                                                              mutation_matrix[pos$Sample,i] <- 1


                                                                          }
                                               else if( length(mutations_gr) <= 0){print(paste0("No intersecting mutations in ", i,"th enhancer"))}

                                           }

return(mutation_matrix)

}

make_CDS_CN_matrix <- function(CN_gr, Gene, CDS_gr,padding=0,normalized=TRUE){
    require(GenomicRanges)
    gene_val <- length(grep(paste0("^",Gene,"$"),CDS_gr$Gene))
    if(gene_val > 0){
        interval_gr <- CDS_gr[grep(paste0("^",Gene,"$"), CDS_gr$Gene)]+padding
        names(interval_gr) <- paste0("exon_",1:length(interval_gr),"_CDS_CN")
        seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
        print(interval_gr)


        CN_matrix <- matrix(0,ncol=length(interval_gr) ,nrow=length(unique(CN_gr$Sample_ID)))

        colnames(CN_matrix) <- names(interval_gr)
        rownames(CN_matrix) <- unique(CN_gr$Sample_ID)
        for(i in colnames(CN_matrix)){
            CN_intersect <- intersect_with_metadata(CN_gr,interval_gr[i])
            if(normalized==TRUE){
                for(j in unique(CN_intersect$Sample_ID)){
                    Sample_mean_CN <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                    Sample_sd_CN <- sd(CN_gr[which(CN_gr$Sample_ID == j)]$CN)
                    CN_val <- (mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN)-Sample_mean_CN)/Sample_sd_CN
                    CN_matrix[j,i] <- CN_val
                }}
            else if(normalized==FALSE){
                for(j in unique(CN_intersect$Sample_ID)){
                    CN_val <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN,na.rm=TRUE)
                    CN_matrix[j,i] <- CN_val
                }}
            print(paste0("Finished processing ",i))}



        } else {print("Cannot find this gene in database")}

    return(CN_matrix)

}

make_enhancer_CN_matrix <- function(CN_gr,interval_gr, blacklist=NULL,normalized=TRUE){
    require(GenomicRanges)

    seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
    if(is.null(blacklist) ==FALSE){
        print("Applying blacklist")
        seqlevelsStyle(blacklist) <- seqlevelsStyle(interval_gr)
        CN_gr <- setdiff_with_metadata(CN_gr,blacklist)
    } else{ print("No blacklist")}
    names(interval_gr) <- paste0(names(interval_gr),"_CN")
    print(interval_gr)
    CN_matrix <- matrix(0,ncol=length(interval_gr) ,nrow=length(unique(CN_gr$Sample_ID)))

        colnames(CN_matrix) <-names(interval_gr)
        rownames(CN_matrix) <- unique(CN_gr$Sample_ID)
        for(i in colnames(CN_matrix)){
            CN_intersect <- intersect_with_metadata(CN_gr,interval_gr[i])
            if(normalized == TRUE){
                for(j in unique(CN_intersect$Sample_ID)){
                Sample_mean_CN <- mean(CN_gr[which(CN_gr$Sample_ID == j)]$CN)

                Sample_sd_CN <- sd(CN_gr[which(CN_gr$Sample_ID == j)]$CN)

                CN_val <- (mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN)-Sample_mean_CN)/Sample_sd_CN

                CN_matrix[j,i] <- CN_val
            }}
            else if(normalized == FALSE){
                for(j in unique(CN_intersect$Sample_ID)){
                    CN_val <- mean(CN_intersect[which(CN_intersect$Sample_ID == j)]$CN,na.rm=TRUE)
                    CN_matrix[j,i] <- CN_val
                }}
            print(paste0("Finished processing ",i))
        }

    return(CN_matrix)

}





calc_model_pval <- function(model){
    f <- summary(model)}
make_l1_enhancer_matrix <- function(reg_element_matrix,expression_matrix,Gene){
    reg_df <- as.data.frame(reg_element_matrix)
    reg_sub <- data.frame(reg_df[intersect(rownames(expression_matrix), rownames(reg_df)),])
    rownames(reg_sub) <- intersect(rownames(expression_matrix), rownames(reg_df))
    colnames(reg_sub) <- colnames(reg_df)
    reg_sub[,Gene] <- NA
    for(i in rownames(reg_sub)){ reg_sub[i,Gene] <- expression_matrix[i,Gene]}
    return(reg_sub)}



make_plot_df <- function(lm_input,Gene){
    out_df <- data.frame(lm_input$model[Gene],lm_input$fitted.values)
    colnames(out_df) <- c("Actual","Estimate")
    out_df$ID <- names(lm_input$fitted.values)
    subtype_column <- grep("Subtype|Status",colnames(lm_input$model))
    if(length(subtype_column) == 0){ out_df$Subtype <- "None"}
    else {out_df$Subtype <- lm_input$model[,subtype_column]}
    final_plot <- ggplot(out_df, aes(x=Actual,y=Estimate,label=ID))+geom_point(aes(color=out_df$Subtype))+geom_text(nudge_y=-0.1,aes(color=out_df$Subtype),size=1.5,check_overlap=TRUE)+geom_abline(linetype="dotted",slope=1,intercept=0)
    return(final_plot)

}

model_AUC <- function(matrix,feature_list,Gene){
    feature_list <- intersect(feature_list, colnames(matrix))
    out_df_final <- data.frame(0,0,0,0)
    colnames(out_df_final) <- c("Actual","Estimate","Feature_Count","Features")
    rownames(out_df_final) <- "Empty"

    for(i in 1:length(feature_list)){
        feature_list_sub <- sample(feature_list,i)
        formula_in <- paste0(feature_list_sub,collapse="+")
        final_formula <- paste0(Gene,"~",formula_in)
        lm_out <- glm(as.formula(final_formula),data=matrix)
        features <- paste(feature_list_sub,collapse="-")
        out_df <- data.frame(lm_out$model[Gene],lm_out$fitted.values,i,features)
        colnames(out_df) <- c("Actual","Estimate","Feature_Count","Features")
        rownames(out_df) <- paste0(names(lm_out$fitted.values),"_",i,"_features")
        out_df_final <- rbind(out_df_final,out_df)
    }
    out_df_final <- out_df_final[2:dim(out_df_final)[1],]
    out_df_final$Status <- NA
    subtype_col <- grep("Status|Subtype",colnames(matrix))
    for(i in rownames(matrix)){index <- grep(paste0(i,"_[0-9]*"),rownames(out_df_final))
                               out_df_final[index,"Status"] <- matrix[i,subtype_col]}

    return(out_df_final)}


feature_selection_experiment <- function(featurelist,matrix,Gene){
    feature_df_list <- list()
    featurelist <- intersect(featurelist,colnames(matrix))
    for (i in c(seq(1,length(featurelist),10),length(featurelist))){
        sample_size <- i
        readout_vec <- c()
        test_features_vec <- c()
        for(j in 1:100){ test_features <- sample(featurelist,sample_size)
                         model <- multiple_regression(matrix,test_features,Gene)
                         readout <- cor(model$model[Gene],model$fitted.values,use='complete.obs')^2
                         readout_vec <- c(readout_vec,readout)
                         test_features_vec <- c(test_features_vec,paste0(test_features,collapse="-"))}
        assign(paste0(i,"_feature_run"), data.frame(readout_vec,test_features_vec))
        feature_df_list <- append(feature_df_list,list(get(paste0(i,"_feature_run"))))
        print(paste0("Finished Processing run with ",i," Features out of a possible ",length(featurelist)))

    }
return(feature_df_list)
}

Collapse_feature_selection <- function(feature_selection_output,feature_description){
    out_df <- data.frame()
    for (i in 1:length(feature_selection_output)){
        sub_df <- feature_selection_output[[i]]
        sub_df$Feature_Count <- length(unlist(strsplit(as.character(sub_df[1,2]),"-")))
        sub_df$Description <- feature_description
        out_df <- rbind(out_df,sub_df)}
    colnames(out_df) <- c("R_Squared","Features","Feature_Count","Description")
    return(out_df)}



extract_best_features <- function(feature_selection_output, feature_description){
    feature_df <- data.frame()
    for(i in 1:length(feature_selection_output)){
        sub <- feature_selection_output[[i]]
        index <- which(sub[,1] == max(sub[,1],na.rm=TRUE))
        features <- as.character(sub[index,2])
        features <- paste0(features,collapse=",")
        out_df <- data.frame(features,feature_description)
        feature_df <- rbind(feature_df,out_df)}
    return(feature_df)}

Group_cross_validation <- function(best_features_output,matrix,num_folds=10,Gene,feature_info){
    in_df <- best_features_output
    final_df <- data.frame()
    for(i in 1:length(in_df[,1])){
        featurelist <- unique(unlist(strsplit(as.character(in_df[i,1]),"-|,")))
        cross_val <- Cross_validation(matrix,num_folds,featurelist,Gene)
        condensed <- Condense_cross_validation_output(cross_val,feature_info)
        final_df <- rbind(final_df,condensed)}
    return(final_df)}


glmnet_regression <- function(featurelist, matrix, Gene, family="gaussian",alpha){
    require(glmnet)
    featurelist_sub <- intersect(featurelist,colnames(matrix))
    check_len <- length(featurelist)-length(featurelist_sub)
    if(check_len > 0){ print("Some of your provided features are *not* in the provided matrix")}
    featurelist <- featurelist_sub
    predictor <- as.matrix(matrix[,featurelist])

    response <- matrix[,Gene]

    fit <- glmnet(x=predictor,y=response,family=family,alpha=alpha)
    return(fit)}

cv_glmnet_regression <- function(featurelist, matrix, Gene, family="gaussian",alpha,nfolds=10, measure='mse'){
    featurelist_sub <- intersect(featurelist,colnames(matrix))
    check_len <- length(featurelist)-length(featurelist_sub)
    if(check_len > 0){ print("Some of your provided features are *not* in the provided matrix")}
    featurelist <- featurelist_sub
    predictor <- as.matrix(matrix[,featurelist])
    response <- matrix[,Gene]
    fit <- cv.glmnet(x=predictor,y=response,family=family,alpha=alpha,type.measure=measure,nfolds=nfolds)
    return(fit)}
run_best_glmnet <- function(cv_glmnet_object, matrix, Gene){

    cols <- coef(cv_glmnet_object)@Dimnames[[1]]
    sub_cols <- intersect(cols, colnames(matrix))
    best_lambda <- cv_glmnet_object$lambda.min
    
    index <- which(cv_glmnet_object$lambda == best_lambda)
    
    num_pos <- cv_glmnet_object$nzero[index]
    newx <- as.matrix(matrix[,sub_cols])

    fit <- predict(cv_glmnet_object,newx=newx, s= best_lambda)
    
    readout <- round(cor(matrix[,Gene], fit)^2,3)
    str(readout)

    return(list(cv_glmnet_object,coef(cv_glmnet_object),readout,num_pos))}


plot_glmnet_regression <- function(glmnet_object,cross_validation_status=FALSE, outfile, title){
    if(cross_validation_status == FALSE){
        pdf(outfile)
        plot(glmnet_object, xvar='lambda')
        plot(glmnet_object, xvar='dev')
        dev.off()
    }
    if(cross_validation_status == TRUE){
        pdf(outfile)
        best_lambda <- glmnet_object$lambda.min
        best_index <- which(glmnet_object$glmnet.fit$lambda == glmnet_object$lambda.min)
        best_model <- glmnet_object$glmnet.fit$dev.ratio[best_index]

        plot(glmnet_object, ylim=c(0,10));abline(v= log(best_lambda),col="red", lty=2,main= title)
        plot(glmnet_object$glmnet.fit, xvar='lambda',label=TRUE);abline(v= log(best_lambda),col="red", lty=2,main= title)
        plot(glmnet_object$glmnet.fit, xvar='dev', xlim= c(0,1),label=TRUE); abline(v= best_model,col="red", lty=2,main= title)
        dev.off()
    }
    return(best_model)
}

plot_case_study <-function(matrix,Gene,title,excluded_columns){
    muts <-matrix[setdiff(colnames(matrix),union(Gene,excluded_columns))]
    mut_rows <- c()
    for(i in rownames(muts)){ mut_row <- muts[i,]; count_element <- as.integer(mut_row !=0); mut_rows <- c(mut_rows, sum(count_element))}
    matrix$Num_mutations <- as.character(mut_rows)
    mut_status <- as.integer(mut_rows !=0)
    mut_status <- gsub(0,"No",mut_status)
    mut_status <- gsub(1,"Yes",mut_status)

    matrix$Mutation_Status <- mut_status

    matrix$Subtype <- expression_df[rownames(matrix),"Status"]


    color_change <- list("LumA"="Blue","LumB"="Green","Basal"="Red","Her2"="Purple","Normal"="White","Unknown"="Black")
    print(paste0("Plotting ", Gene))


    matrix_clean <- muts
    test <- names(table(unlist(unique(apply(matrix_clean,2,unique)))))
    print(test)

    check_cols <- names(table(unlist(unique(apply(matrix_clean,2,unique)))))
    color_pal <- list("-1"="Red","0"="White","1"="Green")


    #colfunc2 <- colorRamp2(breaks= sort(as.integer(check_cols)),colors = unlist(color_pal[check_cols]))




                                        #input_list <- sort(unique(matrix$Num_mutations))
    input_list <- names(which((table(matrix$Num_mutations) >1) == TRUE))

    comparison_list <- list(); for(i in 2:length(input_list)){ int <- list(c("0",input_list[i])); comparison_list <- append(comparison_list,int)}
    heights <- max(matrix[,Gene]) + seq(length.out= length(comparison_list))

                                        #matrix$Num_mutations <- mut_rows
    matrix[is.na(matrix$Subtype),"Subtype"] <- "Unknown"
    if((sum(abs(as.matrix(matrix_clean))) >0) && (dim(as.matrix(matrix_clean))[2] >1)){
        #colfunc2 <- colorRampPalette(colors= unlist(color_pal[check_cols]),bias=10.01, space='rgb')
        #heatmap.2(as.matrix(matrix_clean), col=colfunc2(length(unlist(color_pal[check_cols]))+10),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
                                        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
 #       heatmap.2(as.matrix(matrix_clean), col=colfunc2,RowSideColors=unlist(color_change[matrix$Subtype]), tracecol='white',main=title,cexRow=0.25,cexCol=1.25)


       # heatmap.2(as.matrix(matrix_clean), col=colfunc2,RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)

        p <- ggplot(matrix,aes(x=Num_mutations,y=get(Gene),color=Num_mutations))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Number of Mutated elements",y=paste0(Gene," Expression"))+guides(color=FALSE); p <- p+geom_signif(test='t.test',comparisons=comparison_list,map_signif_level=FALSE,step_increase=0.1,tip_length=0,margin_top=0.4); print(p)
                                                                                                                                                                                          if(length(table(matrix$Num_mutations)) > 2){
            q <- ggplot(matrix,aes(x=Mutation_Status,y=get(Gene),color=Mutation_Status))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Has a mutated element",y=paste0(Gene," Expression")); print(q)} else{print("Not showing T/F for this gene")}
    } else {print("There are no mutations in the gene of interest")}


}

plot_case_study_redux <- function(matrix,Gene,Title,excluded_columns){
    muts <-matrix[setdiff(colnames(matrix),union(Gene,excluded_columns))]
    featurelist <- setdiff(colnames(matrix),excluded_columns)
    lm_out <- multiple_regression(matrix,featurelist,Gene)
    coeffs <- round(summary(lm_out)$coefficients[,1],digits=4)
    coeff_pvals <- round(summary(lm_out)$coefficients[,4],digits=4)

    print(coeffs)
    print(coeff_pvals)

    print(length(abs(coeffs)>0))
    coeff_names <- names(coeffs[which(abs(coeffs) >0)])
    submat <- muts[union(setdiff(coeff_names,"(Intercept)"),grep("_CDS_CN",colnames(matrix),value=TRUE))]
    mut_rows <- c()
    for(i in rownames(muts)){ mut_row <- muts[i,]; count_element <- as.integer(mut_row !=0); mut_rows <- c(mut_rows, sum(count_element))}
    matrix$Num_mutations <- as.character(mut_rows)
    mut_status <- as.integer(mut_rows !=0)
    mut_status <- gsub(0,"No",mut_status)
    mut_status <- gsub(1,"Yes",mut_status)

    matrix$Mutation_Status <- mut_status

    matrix$Subtype <- expression_df[rownames(matrix),"Status"]
    color_change <- list("LumA"="Blue","LumB"="Green","Basal"="Red","Her2"="Purple","Normal"="White","Unknown"="Black")
    print(paste0("Plotting ", Gene))

    matrix_clean <- muts
    check_cols <- names(table(unlist(applyhe(matrix_clean,2,unique))))
    print(check_cols)
    check_cols[which(as.numeric(check_cols) >=2)] <- "2"
    check_cols[which(as.numeric(check_cols) <=-2)] <- "-2"
    check_cols <- unique(check_cols)
    print("New")
    print(check_cols)
    color_pal <- list("-2"="Black","-1"="Red","0"="White","1"="Green","2"="Blue")

   if(length(check_cols) <=4){
                               breaks =seq(as.numeric(check_cols[1]),as.numeric(check_cols[length(check_cols)]),length.out=length(check_cols)+1)} else if(length(check_cols) ==5){ breaks=c(-2,-1,-0.5,0.5,1,2)}

#    breaks <-c(-2,-1,-0.5,0.5,1,2)
    print(breaks)
    print(unlist(color_pal[check_cols]))
    num_breaks=5

    expression_levels <- cut(expression_df[rownames(matrix),Gene],breaks=num_breaks)
    col_levels <- brewer.pal(num_breaks,"RdBu")
    col_list <- list();for(i in names(table(expression_levels))){col_list[[i]] <- col_levels[which(names(table(expression_levels))==i)]}
    rowCols <- unlist(col_list[expression_levels])







                                        #input_list <- sort(unique(matrix$Num_mutations))
    input_list <- names(which((table(matrix$Num_mutations) >1) == TRUE))

    comparison_list <- list(); for(i in 2:length(input_list)){ int <- list(c("0",input_list[i])); comparison_list <- append(comparison_list,int)}
    heights <- max(matrix[,Gene]) + seq(length.out= length(comparison_list))

                                        #matrix$Num_mutations <- mut_rows
    matrix[is.na(matrix$Subtype),"Subtype"] <- "Unknown"
    if((sum(abs(as.matrix(matrix_clean))) >0) && (dim(as.matrix(matrix_clean))[2] >1)){
       #colfunc2 <- colorRampPalette(colors= unlist(color_pal[check_cols]),bias=10.01, space='rgb')
        #heatmap.2(as.matrix(matrix_clean), col=colfunc2(length(unlist(color_pal[check_cols]))+10),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)
                                        #heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]), trace='none',main=title,cexRow=0.25,cexCol=1.25)


        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]),trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256), keysize=1, key.title=NA, key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(matrix_clean), col=unlist(color_pal[check_cols]),RowSideColors=rowCols,trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(submat), col=unlist(color_pal[check_cols]),RowSideColors=unlist(color_change[matrix$Subtype]),trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')
        heatmap.2(as.matrix(submat), col=unlist(color_pal[check_cols]),RowSideColors=rowCols,trace='none', dendrogram='none',main=Title,cexRow=0.05,cexCol=0.65, margins=c(8,4),symkey=FALSE,symbreaks=FALSE,breaks=breaks,sepwidth=c(0.1,0.1),colsep=seq(1,dim(matrix_clean)[2],1),rowsep=seq(1,dim(matrix_clean)[1],1),sepcolor=rgb(col2rgb('gray91')[1,],col2rgb('gray91')[2,],col2rgb('gray91')[3,],alpha=232,maxColorValue=256),keysize=1,key.title=NA,key.xlab='Mutation Score',density.info='none')


        p <- ggplot(matrix,aes(x=Num_mutations,y=get(Gene),color=Num_mutations))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Number of Mutated elements",y=paste0(Gene," Expression"))+guides(color=FALSE); p <- p+geom_signif(test='t.test',comparisons=comparison_list,map_signif_level=FALSE,step_increase=0.1,tip_length=0,margin_top=0.4); print(p)
      if(length(table(matrix$Num_mutations)) > 2){
            q <- ggplot(matrix,aes(x=Mutation_Status,y=get(Gene)))+geom_boxplot(width=0.5)+geom_jitter(color=unlist(color_change[matrix$Subtype]),width=0.2)+labs(x="Has a mutated element",y=paste0(Gene," Expression"))+geom_signif(test='t.test', comparisons=list(c("No","Yes")),map_signif_level=FALSE,tip_length=0); print(q)} else{print("Not showing T/F for this gene")}}             else {print("There are no mutations in the gene of interest")}


}



multiple_regression <- function(matrix,feature_list,Gene){
    feature_list <- intersect(feature_list, colnames(matrix))

    formula_in <- paste0(feature_list,collapse="+")
    str(formula_in)
    final_formula <- paste0(Gene,"~",formula_in)
    lm_out <- bglm(as.formula(final_formula), data=matrix,maxit=100)
    return(lm_out)}

make_conting_table_model_AUC <- function(model_AUC_output,num_steps=20){
    conting_table_final <- data.frame(0,0,0,0,0,0,0,"Empty")
    colnames(conting_table_final) <- c("TP","FP","TN","FN","Feature_Count","TPR","FPR","Features")
    for(i  in unique(model_AUC_output$Feature_Count)){
        in_df <- model_AUC_output[which(model_AUC_output$Feature_Count == i),]
        step_counter <- seq.int(min(in_df$Estimate),max(in_df$Estimate),length.out=num_steps)

        TP_vec <- c()
        FP_vec <- c()
        TN_vec <- c()
        FN_vec <- c()

        readout_list <- list("11"="TP","10"="FP","01"="FN","00"="TN")

        for (j in step_counter){
            Prediction <- as.integer(in_df$Estimate > j)
            Actual <- as.integer(in_df$Status %in% c("LumA","LumB"))
            Comparison <- as.integer(Prediction == Actual)




            readout <- paste0(Prediction,Actual)
            readout_values <- table(unlist(readout_list[readout]))
            TP_vec <- as.vector(c(TP_vec,readout_values["TP"]))
            FP_vec <- as.vector(c(FP_vec,readout_values["FP"]))
            TN_vec <- as.vector(c(TN_vec,readout_values["TN"]))
            FN_vec <- as.vector(c(FN_vec,readout_values["FN"]))

            }
        conting_table <- data.frame(TP_vec,FP_vec,TN_vec,FN_vec,i)
        colnames(conting_table) <- c("TP","FP","TN","FN","Feature_Count")
        conting_table[is.na(conting_table)] = 0
        conting_table$TPR <- conting_table$TP/(conting_table$TP+conting_table$FN)
        conting_table$FPR <- 1 -(conting_table$TN/(conting_table$TN+conting_table$FP))
        conting_table$Features <-  as.character(unique(in_df$Features))



        conting_table_final <- rbind(conting_table_final,conting_table)

    }
    conting_table_final <- conting_table_final[2:dim(conting_table_final)[1],]
    conting_table_final$Features <- as.character(conting_table_final$Features)

    return(conting_table_final)}

plot_conting_table <- function(conting_table_AUC, outfile){
    pdf(outfile)
    if(length(unique(conting_table_AUC$Feature_Count)) > 5 ){
        subsample <- sample(unique(conting_table_AUC$Feature_Count),5)}
    else if(length(unique(conting_table_AUC$Feature_Count)) < 5){
        subsample <- sample(unique(conting_table_AUC$Feature_Count),length(unique(conting_table_AUC$Feature_Count)))}

    conting_table_sub <- conting_table_AUC[which(conting_table_AUC$Feature_Count %in% subsample),]
    AUC_vec <- c()
    for(i in conting_table_sub$Feature_Count){sub <- conting_table_sub[which(conting_table_sub$Feature_Count == i),];     AUC <- integrate.xy(sub$FPR,sub$TPR); AUC_vec <- c(AUC_vec,AUC)}
    AUC_vec <- round(AUC_vec,digits=3);str(AUC_vec);str(conting_table_sub)
    conting_table_sub$Feature_Count <- paste0(conting_table_sub$Feature_Count,":AUC=",AUC_vec)
    p <- ggplot(conting_table_sub,aes(x=FPR,y=TPR,color=Feature_Count))+geom_line()+geom_abline(slope=1,intercept=0,linetype="dotted")
    q <- p+labs(color="Feature Count",caption= gsub(".pdf","",outfile))
    print(q)

    dev.off()
    out <- gsub(".pdf",".txt",outfile)
    write.table(conting_table_sub,out,sep='\t',quote=F,row.names=F)
    return(conting_table_sub)
}

Plot_vector_as_bar_chart <- function(vector,outfile,x_field,y_field,negative_axis=FALSE,return_df=FALSE){
    df <- data.frame(vector)
    df$Name <- names(vector)
    colnames(df) <- c("Value","Name")

    str(df)
    vector_mean <- mean(abs(df$Value))
    print(vector_mean)
    if(negative_axis ==TRUE){
        p <- ggplot(df,aes(x=Name,y=Value,label=Name))+labs(x=x_field,y=y_field)+geom_col(alpha=0.65)+geom_text(nudge_y=0.01,size=3,check_overlap=TRUE)+theme(axis.text.x=element_text(face='bold',size=3,angle=45))+geom_hline(yintercept=vector_mean,linetype='dotted',col="salmon",size=2)+geom_hline(yintercept=-1*vector_mean,linetype='dotted',col="salmon",size=2)}
    else if(negative_axis==FALSE){
        p <- ggplot(df,aes(x=Name,y=Value,label=Name))+labs(x=x_field,y=y_field)+geom_col(alpha=0.65)+geom_text(nudge_y=0.01,size=3,check_overlap=TRUE)+theme(axis.text.x=element_text(face='bold',size=3,angle=45))+geom_hline(yintercept=vector_mean,linetype='dotted',col="salmon",size=2)}
    pdf(outfile)
    print(p)
    dev.off()
    if(return_df == TRUE){ return(df)} else{print("Finished")}
}


prepare_matrix <- function(matrix, gene) {
    pheno_input <- matrix[,gene]
    cleaned_matrix <- matrix
    cleaned_matrix[,gene] <- NULL
    return(list(pheno_input,cleaned_matrix))
}

pval_calculator <- function(pheno_input, x_input ){


    X_mx <- as.matrix(cbind(1,x_input))
    n_samples <- length(x_input)[1]


    MLE_beta <- ginv(t(X_mx) %*% X_mx) %*% t(X_mx) %*% pheno_input

    y_hat <- X_mx %*% MLE_beta


    SSM <- sum((y_hat - mean(pheno_input))^2)
    SSE <- sum((pheno_input - y_hat)^2)

    df_M <- 2
    df_E <- n_samples - 3
    MSM <- SSM / df_M
    MSE <- SSE / df_E


    Fstatistic <- MSM / MSE
    pval <- pf(Fstatistic,df1 = 2, n_samples-3,lower.tail = FALSE)



    return(pval)
}
significant_pval_ID <- function(pval_df){ cutoff <- 0.05/length(pval_df[,1]); significant_ids <- which(pval_df[,1] <= cutoff); significant <- pval_df[significant_ids,]; names(significant) <- rownames(pval_df)[significant_ids]; return(significant)}

table_to_granges <- function(table,strand=NULL){ require(GenomicRanges)
                                                 if(is.null(strand)){
                                                     gr <- GRanges(seqnames=Rle(table[,1]), ranges= IRanges(start=as.numeric(table[,2]), end=as.numeric(table[,3])),strand="*")
                                                     ##print(setdiff(colnames(table),x=colnames(table)[4:ncol(table)]))
                                                     if(ncol(table) >=4){
                                                     gr@elementMetadata@listData <- as.list(table[setdiff(colnames(table)[1:3],x=colnames(table)[4:ncol(table)])])} 
                                                 } else {gr <- GRanges(seqnames=Rle(table[,1]), ranges= IRanges(start=table[,2], end=table[,3]),strand=table[,strand])
                                                                                                                                                               gr@elementMetadata@listData <- as.list(table[setdiff(grep("strand", colnames(table),ignore.case=T),x=4:ncol(table))])}
                        return(gr)}

tidy_matrix <- function (matrix){
    matrix <- as.data.frame(lapply(matrix,function(y) as.numeric(gsub(-Inf,-10000,y))))
    matrix[is.na(matrix)] <- 0
    return(matrix)

}

violin_plots <- function(matrix, list_of_features, group_column=NULL,outfile){
    pdf(outfile)
    featurelist <- sort(intersect(list_of_features,colnames(matrix)))
    str(featurelist)
    str(list_of_features)
    print(setdiff(list_of_features,featurelist))
    if(is.null(group_column) == FALSE){group_column <- intersect(group_column, colnames(matrix))}
    else {group_column <- group_column}

    submatrix <- matrix[,c(featurelist,group_column)]
    for (i in featurelist){ p <- ggplot(submatrix, aes(x=get(group_column),y=get(i), fill=get(group_column)))+geom_violin(trim=FALSE)+labs(fill=group_column)+xlab(group_column)+ylab(i)
                            q <- p+stat_summary(fun.data="mean_sdl", geom="crossbar",width =0.04)+geom_signif(test='t.test',comparisons = list(c("LumA","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+1.5)+geom_signif(test='t.test',comparisons = list(c("LumB","Normal")), map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+1) + geom_signif(test='t.test',comparisons= list(c("Her2","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2) + geom_signif(test='t.test',comparisons = list(c("Basal","Normal")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2.5)+ geom_signif(test='t.test',comparisons = list(c("Yes","No")),map_signif_level=TRUE,hjust=1,y_position = max(submatrix[,i])+2.5)
                            print(q)
                            print(i)}
    dev.off()}



get_CN_expression_correlation <- function(genelist,matrix){
    cor_vec <- c()
    for(i in genelist){ expression <- matrix[i]
                        CN <- get(paste0("CNA_table_",i))$Copy_Number
                        cor_val <- cor(expression,CN,use='complete.obs')
                        cor_vec <- c(cor_vec,cor_val)}
    names(cor_vec) <- genelist
    return(cor_vec)}


Copy_Number_Eval_all <- function(genelist, expression_matrix, Gene, CN_gr, CDS_gr,padding=0 ){
    require(GenomicRanges)
    CNA_table_final <- as.data.frame(expression_matrix[,Gene])
    colnames(CNA_table_final) <- Gene
    for(i in genelist){
        gene_val <- length(grep(paste0("^",i,"$"),CDS_gr$Gene))
        if(gene_val > 0){ interval_gr <- CDS_gr[grep(paste0("^",i,"$"), CDS_gr$Gene)[1]]+padding
                          samples_intersecting_Gene <- intersect_with_metadata(CN_gr,interval_gr)
                          seqlevelsStyle(CN_gr) <- seqlevelsStyle(interval_gr)
                          rownames(CNA_table) <- rownames(expression_matrix)
                          CNA_table$Copy_Number <- NA
                          colnames(CNA_table) <- c(Gene,"Copy_Number")
                          CNA_table$ID <-  rownames(CNA_table)
                          for(j in rownames(CNA_table)){ index <- which(samples_intersecting_Gene$Sample_ID == j)
                                                         CN <- mean(samples_intersecting_Gene$CN[index])
                                                         CNA_table[j,"Copy_Number"] <- CN }
                          CNA_table$Subtype <- NA
                          for(k in rownames(CNA_table)){
                              sample_subtype <- expression_matrix[k,"Status"]
                              CNA_table[k,]$Subtype <- sample_subtype}
                          colnames(CNA_table) <- c(Gene,paste0(i,"_Copy_Number"),"ID","Subtype")
                          CNA_table_final <- cbind(CNA_table_final, CNA_table[,paste0(i,"_Copy_Number")])
                          col_names <- colnames(CNA_table_final)
                          index <- grep("CNA_table",col_names)
                          col_names[index] <- paste0(i,"_CN")
                          colnames(CNA_table_final) <- col_names
                          print(paste0("Finished ",i))

                      }
        else if(gene_val ==0){
                                        #print(paste0("No CDS for ",i))
        }

    }
    CNA_table_final$ID <-CNA_table$ID
    CNA_table_final$Subtype <- CNA_table$Subtype
return(CNA_table_final)
}


get_regulator_TFs <- function(graph,Gene){
    require(igraph)
    index <- which(names(V(graph)) == Gene)
    if(length(index) >0 ){    regulators <- names(which(graph[,index]!=0))} else{ regulators <- ""}
    return(regulators)}

get_regulatory_targets <- function(graph,TF){
    require(igraph)
    index <- which(names(V(graph)) == TF)
    if(length(index) > 0){    targets <- names(which(graph[index,]!=0))} else {targets <- ""}
    return(targets)}

subset_graph <- function(graph, subset,mode=NULL,TF_filter=NULL){
    require(igraph)
    require(dplyr)

    el <- as.data.frame(as_edgelist(graph))


    if(is.null(mode)==TRUE){
        if(is.null(TF_filter) ==TRUE){
           sub_el <- rbind(filter(el, V1 %in% subset),filter(el, V2 %in% subset))
           out <- graph_from_edgelist(as.matrix(sub_el), directed=TRUE)} else if (is.null(TF_filter)==FALSE){
               print("double_subset")
               print(subset)
               sub_el <- rbind(filter(el, V1 %in% intersect(subset,unique(el[,1]))),filter(el, V2 %in% intersect(subset,unique(el[,1]))))
               sub_el <- filter(sub_el, V2 %in% unique(el[,1]))
               out <- graph_from_edgelist(as.matrix(sub_el), directed=TRUE)
           }

   } else if (mode=="out"){
       if(is.null(TF_filter) ==TRUE){
           sub_el <- rbind(filter(el, V1 %in% subset),filter(el, V2 %in% subset))
           out <- graph_from_edgelist(as.matrix(sub_el), directed=TRUE)} else if (is.null(TF_filter)==FALSE){
               
               sub_el <- rbind(filter(el, V1 %in% intersect(unique(el[,1]),subset)),filter(el, V2 %in% unique(el[,1])))
               out <- graph_from_edgelist(as.matrix(sub_el), directed=TRUE)
                           ##sub_el <- as.matrix(filter(sub_el, V1 %in% subset))
            ##out <- graph_from_edgelist(sub_el, directed=TRUE)}
           }} else if(mode=="in"){
                sub_el <- as.matrix(filter(el, V2 %in% subset))
                out <- graph_from_edgelist(sub_el, directed=TRUE)}
    
    return(out)}

delete.na <- function(DF, n=0) {
    DF[rowSums(is.na(DF)) <= n,]
}
split_vector <- function(vector, split_size){ split(vector, ceiling(seq_along(vector)/split_size))}

RIGHT <- function(x,n){
    substring(x,nchar(x)-n+1)
}
get_hub_TFS <- function(graph,q=0.75){
    outdegree <- igraph::degree(graph,mode='out')
    outdegree <- outdegree[which(outdegree >0)]
    hubs <- sort(names(outdegree[which(outdegree >= quantile(outdegree,q))]))
    return(hubs)}

get_common_hub_TFs <- function(hublist, graphs_to_exclude) {
    sub_list <- hublist[setdiff(names(hublist),graphs_to_exclude)]
    common_vec <- sub_list[[1]]
    for(i in 2:length(sub_list)){
        common_vec <- intersect(common_vec,sub_list[[i]])}
    blacklist <- unique(unlist(hublist[graphs_to_exclude]))
    common_vec <- setdiff(common_vec,blacklist)
    return(common_vec)}


importJaspar <- function(file=myloc) {
    vec <- readLines(file)
    vec <- gsub("\t"," ",vec)
    vec <- gsub("\\[|\\]", "", vec)
    start <- grep(">", vec); end <- grep(">", vec) - 1
    pos <- data.frame(start=start, end=c(end[-1], length(vec)))
    pwm <- sapply(seq(along=pos[,1]), function(x) vec[pos[x,1]:pos[x,2]])
    pwm <- sapply(seq(along=pwm), function(x) strsplit(pwm[[x]], " {1,}"))
    pwm <- lapply(seq(along=start), function(x) matrix(as.numeric(t(as.data.frame(pwm[(pos[x,1]+1):pos[x,2]]))[,-1]), nrow=4, dimnames=list(c("A", "C", "G", "T"), NULL)))
    names(pwm) <- gsub(">", "", vec[start])
    return(pwm)
}

Jaspar2meme <- function(encode_motifs_file, outdir){
    in_motifs <- importJaspar(encode_motifs_file)
    out_motifs <- list(); for(i in 1:length(in_motifs)){ motif <- t(in_motifs[[i]])
                                                         motif <- motif/rowSums(motif)
                                                         motif_name <- gsub(" ","~",names(in_motifs[i]))
                                                         outlist <- list(motif)
                                                         outname <- paste0(outdir, motif_name,"_MEME.txt")

                                                         writeLines("MEME version 4", outname)
                                                         write("ALPHABET= ACGT",outname,append=T)
                                                         write("strands: +-",outname,append=T)
                                                         write(paste0("MOTIF ", motif_name), outname,append=T )
                                                         write(paste0("letter-probability matrix: alength=4 w= ",dim(motif)[1]),outname,append=T)
                                                         vec <- c()
                                                         for(j in 1:dim(motif)[1]){
                                                             out_vec <- paste(motif[j,],collapse=" ")
                                                             vec <- c(vec,out_vec)}
                                                         write(vec, outname, append=TRUE)}

}

FIMO_to_granges <- function(FIMO_bed){
    require(GenomicRanges)
    bed <- read.table(FIMO_bed, header=T, stringsAsFactors=F,sep='\t')
    gr <- GRanges(seqnames=Rle(bed[,1]), ranges= IRanges(start=bed[,2], end=bed[,3]),strand=bed[,4], PWM_Score=bed[,5],P_value=bed[,6], sequence=bed[,7], motif=bed[,8])
                                 return(gr)}

Collapse_gr_keep_Gene <- function(gr){
    require(GenomicRanges)
    reduced <- reduce(gr)
    final_gr <- reduced
    for (i in 1:length(reduced)){
        overlap_enhancer = findOverlaps(reduced[i], gr)
        row_names = vector()
        for( j in 1:length(overlap_enhancer)){
            names = gr[overlap_enhancer[j]@to]$Gene
            row_names = c(names, row_names)

        }

        final_gr$nearby[i] = row_names
    }
    return(final_gr)}

make_qqplot_df <- function(pval_vector){
    pval_df <- data.frame(sort(-log10(pval_vector)),stringsAsFactors=F)
        colnames(pval_df) <- "observed"

    pval_df$expected <- sort(-log10(runif(length(pval_vector))))
    return(pval_df)
}
qqplot <- function(pval_vector){
        make_qqplot_df(pval_vector)
        p <- ggplot(pval_df,aes(expected,observed))+geom_point(size=2)+geom_abline(slope=1,intercept=0,size=2,color="red",linetype=2, alpha=0.5)+labs(title=paste0("QQplot"))
        return(p)}

lseq <- function(from=1, to=100000, length.out=19){ exp(seq(log(from), log(to), length.out=length.out))}

get_nearest_gr <- function(gr1,gr2){
    require(GenomicRanges)
    int <- distanceToNearest(gr1,gr2)
    sub <- gr1[int@from]
    sub$Distance <- int@elementMetadata@listData$distance
    sub$closest <- int@to
    return(sub)}


pintersect_with_metadata <- function(gr1,gr2){
    require(GenomicRanges)
    hits <- findOverlaps(gr1,gr2)
    gr.over <- pintersect(gr1[queryHits(hits)],gr2[subjectHits(hits)])
    gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x) sum(width(x)))
    gr1$overlap<- 0
    gr1$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
    return(gr1)}

empty_gr <- function(){
    require(GenomicRanges)
    gr <- GRanges(seqnames='chr1',ranges=IRanges(start=1,end=2))
    return(gr)}

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

run_tsne <- function(df,sample_margin="column",perplexity=NULL,...){
    require(Rtsne)
    if(sample_margin == "column"){ df <- as.data.frame(t(df),stringsAsFactors=F)
                                   rownames(df) <- gsub("\\.","-",rownames(df))
                               }
    else if (sample_margin == "row") { df <- df}
    if(is.null(perplexity) == TRUE){ hyperparam = round((nrow(df)-1)/3)-1}
    else{ hyperparam=perplexity}
    out <- Rtsne(df,check_duplicates=F,perplexity=hyperparam,...)$Y
    #print(Rtsne(df,check_duplicates=F,perplexity=hyperparam,...)$theta)
    out <- as.data.frame(out,stringsAsFactors=F)
    out$Sample <- rownames(df)
    rownames(out) <- out[,3]
    colnames(out) <- c("Dim1","Dim2","Sample")
    return(out)
}

get_metadata_tsne <- function(df, metadata_file,column=2,col_name="Subtype"){
    df[col_name] <- metadata_file[rownames(df),column]
    return(df)}
get_clusters_tsne <- function(tsne_df,method=c("hclust","dbscan","kmeans"),cutoff_quantile=0.95,opt_count=1,minPoints=3,eps=NULL,kmin_samples=3,keep_outliers=TRUE,outlier_method=c("Centroid","Agglomerative"),kmean_clusters=NULL){
    require(dplyr)
    sub_df <- tsne_df[,grep("Dim", colnames(tsne_df))]
##    str(sub_df)
    sub_df$Sample <- tsne_df$Sample
    in_df <- sub_df

    if(method == "hclust"){
    clusters_tsne <- hclust(dist(in_df[,c(1,2)]))
    members <- cutree(clusters_tsne, h=quantile(clusters_tsne$height,cutoff_quantile))
 ##   str(members)

        print("Clusters identified")
    cluster_out <- as.data.frame(cbind(members,in_df$Sample),stringsAsFactors=F)
    #str(in_df$Sample)
    #str(cluster_out)
    colnames(cluster_out) <- c("Cluster","Sample")
        tsne_out <- merge(tsne_df,cluster_out,by="Sample")
        print("Merging clusters to tSNE")
        rownames(tsne_out) <- tsne_out$Sample
    colnames(tsne_out) <- gsub("members","Cluster",colnames(tsne_out))

    } else if(method=="kmeans"){

        if(is.null(kmean_clusters)){
        optimal_k <- get_optimal_k_clusters(tsne_df,method="kmeans")
        ##print(optimal_k)
        print(paste0("Optimal K=",optimal_k[[1]]))
        optimal_k[[1]] <- ifelse(optimal_k[[1]]==1,3,optimal_k[[1]])
        print(optimal_k[[1]])
 ##       str(sub_df)
        tsne_df$Cluster <- as.character(kmeans(sub_df[1:2],optimal_k[[1]])$cluster)} else{
                                                                                       tsne_df$Cluster <- as.character(kmeans(sub_df[1:2],kmean_clusters)$cluster)}

   ##     str(tsne_df)
        tsne_out <- tsne_df
        ## stop()


        } else if(method == "dbscan"){
    require(dbscan)
    require(dplyr)
    if(is.null(opt_count)){
        if(is.null(eps)){
            opt_eps <- get_optimal_dbscan_eps(in_df,opt_count=1,min_samples=minPoints,kmin_samples=kmin_samples)} else{opt_eps <- eps}
    } else if(!is.null(opt_count)) { if(is.null(eps)){ opt_eps <- get_optimal_dbscan_eps(in_df[1:2],opt_count=opt_count,min_samples=minPoints,kmin_samples=kmin_samples)} else { opt_eps <- eps}}

    members <- dbscan(in_df[1:2],opt_eps,minPoints)$cluster
    print(table(members))
    if(any(members == 0)){ index <- which(members==0)
                           print("Working on outliers")

                           if(keep_outliers==TRUE){
                               print("Preserving outliers")
                               addition <- 1:length(index)
                           members[index] <- max(members)+addition

                           tsne_df$Cluster <- members} else if (keep_outliers==FALSE){
                               print("Merging outliers")
                               if(outlier_method=="Centroid"){
                               sub_tsne <- tsne_df
                               sub_tsne$Cluster <- members
                               sub_tsne_real <- dplyr::filter(sub_tsne, Cluster !=0)
                               if(nrow(sub_tsne_real) <1){ tsne_df$Cluster <- 0} else{
                               centroid_df <- sub_tsne_real %>% group_by(Cluster) %>% dplyr::select(Dim1,Dim2) %>% summarize_all(mean) %>% ungroup

                               sub_tsne_outlier <- dplyr::filter(sub_tsne, Cluster ==0)
                               sub_tsne_outlier$Cluster <- as.numeric(unlist(centroid_df[apply(sub_tsne_outlier,1, function(x) get_closest(x[1:2], centroid_df[2:3])),"Cluster"]))
                               #str(sub_tsne_outlier)
                               tsne_df <- rbind(sub_tsne_outlier,sub_tsne_real)
                               tsne_df <- tsne_df[match(in_df$Sample,tsne_df$Sample),]
##                               str(tsne_df)
                               tryCatch( {rownames(tsne_df)= tsne_df$Sample} ,error=function(e) {rownames(tsne_df) =1:nrow(tsne_df)})
                           }} else if(outlier_method=="Agglomerative"){

                               sub_tsne <- tsne_df;
                               sub_tsne$Cluster <- members
                               sub_tsne_real <- dplyr::filter(sub_tsne, Cluster !=0)

                               sub_tsne_outlier <- dplyr::filter(sub_tsne, Cluster ==0)
                               print(paste0("Merging any samples within range ", eps))
                               sub_tsne_outlier$Cluster <- dbscan(sub_tsne_outlier[1:2],opt_eps, minPts=2)$cluster

                               final_non_outlier <- dplyr::filter(sub_tsne_outlier, Cluster!=0)
                               final_non_outlier$Cluster <- final_non_outlier$Cluster+max(sub_tsne_real$Cluster)

                               sub_tsne_real <- rbind(sub_tsne_real,final_non_outlier)
                               centroid_df <- sub_tsne_real %>% group_by(Cluster) %>% dplyr::select(Dim1,Dim2) %>% summarize_all(mean) %>% ungroup

                               final_outlier <- dplyr::filter(sub_tsne_outlier, Cluster==0)
                               final_outlier$Cluster <- as.numeric(unlist(centroid_df[apply(final_outlier,1, function(x) get_closest(x[1:2], centroid_df[2:3])),"Cluster"]))
                               tsne_df <- rbind(final_outlier,sub_tsne_real)
                               tsne_df <- tsne_df[match(in_df$Sample,tsne_df$Sample),]
                                tryCatch( {rownames(tsne_df)= tsne_df$Sample} ,error=function(e) {rownames(tsne_df) =1:nrow(tsne_df)})


                            } else { print(paste0("Outlier method:",outlier_method," not implemented. Options are Centroid/Agglomerative"))
                                stop()
                            }
                                                     }}


##    str(tsne_df)
     ## stop()

    tsne_out <- tsne_df

        }
    
    tsne_out$Cluster <- factor(tsne_out$Cluster,levels=sort(as.numeric(sort(unique(tsne_out$Cluster)))))
    print(paste0("Found ",length(unique(tsne_out$Cluster))," clusters"))

        return(tsne_out)}



plot_tsne <- function(df,title,color="Subtype",label="Subtype",viz=c("basic","medium","enhanced"),shape=NULL,alpha=0.5,size=3,font_size=3,legend_col=2){
    require(ggplot2)
#    require(ggforce)
    require(ggrepel)
    require(dplyr)

    if(color != label){
        if(is.null(shape)==TRUE){
            centroid_df <- df %>% group_by(get(color),get(label)) %>% dplyr::select(Dim1,Dim2) %>% summarize_all(mean) %>% ungroup
            colnames(centroid_df)[1:2] <- c(color,label) } else {
                                                             ##centroid_df <- df %>% group_by(get(color),get(label),get(shape)) %>% dplyr::select(Dim1,Dim2) %>% summarize_all(mean) %>% ungroup
                                                               ##   colnames(centroid_df)[setdiff(1:ncol(centroid_df),grep("Dim",colnames(centroid_df)))] <- unique(c(color,label,shape))
                                                                 ## print(str(centroid_df))
                                                              }
    } else { centroid_df <- df %>% group_by(get(label)) %>% dplyr::select(Dim1,Dim2) %>% summarize_all(mean) %>% ungroup
                                                         colnames(centroid_df)[1] <- label
                                                     }
    if(viz=="enhanced"){
        p <- ggplot(df, aes(Dim1,Dim2, fill=!!sym(color),label=!!sym(label)))+geom_point(size=size, alpha=alpha,shape=21,color="gray50")+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=1))+geom_mark_rect(aes(label=!!sym(label),fill=!!sym(color)),show.legend=FALSE,label.buffer=unit(1,"mm"),expand=unit(3,"mm"),label.fontsize=10,con.type="straight")} else if(viz=="basic"){
            p <- ggplot(df, aes(Dim1,Dim2, fill=!!sym(color),label=!!sym(label)))+geom_point(size=size,shape=21,color="gray50",stroke=0.1,alpha=alpha)+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=legend_col))+geom_text(fontface="bold",size=font_size,check_overlap=TRUE,show.legend=FALSE)} else if(viz=="medium"){
                if(is.null(shape)){ p <- ggplot(df, aes(Dim1,Dim2,fill=!!sym(color)))+geom_point(size=size,alpha=alpha,shape=21, color="gray50")+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=legend_col))+geom_text_repel(data=centroid_df, aes(label=!!sym(label)), segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE)}
                else if(!is.null(shape)){p <- ggplot(df, aes(Dim1,Dim2,label=!!sym(label),fill=!!sym(color)))+geom_point(aes(fill=!!sym(color),shape=!!sym(shape)),size=size,alpha=alpha,color="gray50")+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=FALSE,color=FALSE)+scale_shape_manual(values=c(21:25))+geom_text_repel(data=centroid_df, aes(label=!!sym(label),fill=!!sym(color)), segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE)

                    ###+geom_text_repel(data=centroid_df, aes_string(label=label), segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE)
                }
            }




    return(p)
}

get_drivers_from_clusters <- function(df,tsne_df, rank=TRUE,metadata_col="Cluster",subset=TRUE,qval_cutoff=0.2,ratio_cutoff=NULL){
    require(reshape2)
    require(dplyr)
    tsne_df$Cluster_drivers <- "None"
    tsne_df$Signif_drivers <- "None"
    tsne_df$Top_Driver <- "None"
    tsne_df$Top_5 <- "None"
    colnames(df) <- gsub("\\.","-",colnames(df))
    if(subset==TRUE){ df <- df[,which(colnames(df) %in% tsne_df$Sample)]
                      tsne_df <- tsne_df[which(tsne_df$Sample %in% colnames(df)),]
                  } else { print("Not filtering")}

    pval_all <- data.frame(stringsAsFactors=F)
    print("Initializing")
 
    for(i in sort(unique(tsne_df[,metadata_col]))){
        sample_list <- tsne_df[which(tsne_df[,metadata_col] == i),"Sample"]
#        str(sample_list)
        mat_sub <- df[sample_list]
        other_mat <- df[setdiff(colnames(df),sample_list)]
        pval_df <- data.frame(stringsAsFactors=F)
                                        #       str(sample_list)
##        str(mat_sub)
                                        #        str(other_mat)

        
        for(j in rownames(other_mat)){
            if(length(sample_list) >1){
            pval <- wilcox.test(unlist(mat_sub[j,]),unlist(other_mat[j,]))$p.value} else if(length(sample_list)==1){ pval <- 0.05}

            cluster_avg <- mean(unlist(mat_sub[j,]))
            alt_avg <- mean(unlist(other_mat[j,]))
            sample_num <- length(unlist(mat_sub[j,]))
            cluster_id <- i
            ratio <- mean(unlist(mat_sub[j,]))/mean(unlist(other_mat[j,]))

#            str(unlist(ratio))
            if(rank == TRUE){ ratio <- 1/ratio
                sign <- ratio > 1 } else if(rank ==FALSE){ratio <- ratio
                                                  sign <- ratio > 1 & cluster_avg >= alt_avg}



            pval_df <- rbind(pval_df,cbind(cluster_avg,alt_avg,sample_num,cluster_id,pval,j,sign,ratio),stringsAsFactors=F)}
                                        #str(pval_df)
        print(paste0("Finished cluster ",i))
#        str(pval_df)
        pval_df[,c(1:5,8)] <- apply(pval_df[,c(1:5,8)],2, function(x) as.numeric(x))
        pval_df$qval <- p.adjust(pval_df$pval,method="BH")
        pval_df[is.na(pval_df)] <- 0

        if(all(pval_df$qval >=qval_cutoff) == TRUE){
            pval_df$qval <- pval_df$pval
            pval_df$Flag <- "FDR_Fail"} else {
                pval_df$Flag <- " "}

#        print(table(pval_df$sign))
#        pval_df <- pval_df[which(pval_df$sign =="TRUE"),]
        pval_df <- pval_df[order(-as.numeric(as.logical(pval_df$sign)),pval_df$qval,-pval_df$ratio),]

        ##Modified to use ratio as well as sign of changes##

#        top10_df <- pval_df[1:10,]
#        top10 <- top10_df[,"j"]

########################################################################

        if(is.null(ratio_cutoff)){
            signif_df <- pval_df[which(pval_df$qval <= qval_cutoff & pval_df$sign == "TRUE" & pval_df$ratio >= quantile(pval_df$ratio,0.95)),]}
        else{signif_df <- pval_df[which(pval_df$qval <= qval_cutoff & pval_df$sign == "TRUE" & pval_df$ratio >= ratio_cutoff),]}
        signif_df <- signif_df[order(-signif_df$ratio),]
#        str(signif_df)

        top10_df <- signif_df[1:10,]
        top10 <- top10_df[,"j"]
        top10 <- gsub("-NA","",top10)
        top5 <- top10[1:5]
        ##        str(top10_df)
##        print(top10)
        print(top5)

        top_driver <- top10_df[which(top10_df$ratio == max(top10_df$ratio,na.rm=TRUE)),"j"]
##        str(top_driver)

#        str(tsne_df[sample_list,])
#        str(rownames(tsne_df))
        tsne_df[sample_list,"Cluster_drivers"] <- paste0(top10,collapse="-")
        tsne_df[sample_list,"Signif_drivers"] <- paste0(signif_df[,"j"],collapse="-")
        tsne_df[sample_list,"Top_Driver"] <- paste0(top_driver,collapse="-")
        tsne_df[sample_list,"Top_5"] <- paste0(top5,collapse="-")





        pval_df$Cluster <- paste0("Cluster-",i)
        pval_all <- rbind(pval_all,pval_df,stringsAsFactors=F)
       # str(pval_all)
    }
    tsne_df$Plot <- paste0("Clust-",tsne_df$Cluster,":",tsne_df$Top_5)
    tsne_df <- tsne_df %>% group_by(Cluster) %>% mutate(Disease=paste0(paste0(unique(Subtype),collapse="-"),"~Clust-",Cluster)) %>% data.frame
    tsne_df[c("Cluster_drivers","Signif_drivers")] <- apply(tsne_df[c("Cluster_drivers","Signif_drivers")],2, function(x) gsub("-NA","",x))
##    str(tsne_df)
##    str(tsne_df$Sample)
    rownames(tsne_df) <- tsne_df$Sample
    return(list(pval_all,tsne_df))
}

make_clusters_from_metric_list <- function(metric_list,cutoff_quantile=0.95,metadata_file,perp=10,method="dbscan"){
    metric_clusters <- list()
    for(i in 1:length(metric_list)){
        print(paste0("Working ",names(metric_list)[i]))
        df <- metric_list[[i]]

        tsne_out <- run_tsne(df,perplexity=perps)
        print("tSNE done")

        tsne_out <- get_clusters_tsne(tsne_out,method)

        tsne_out <- get_metadata_tsne(tsne_out, metadata_file)
        str(rownames(tsne_out))

        #print(head(tsne_out))
        tsne_out <- get_drivers_from_clusters(df,tsne_out)
        #str(tsne_out)
        metric_clusters <- c(metric_clusters, list(tsne_out)); print(paste0("Finished metric ",i))}
    names(metric_clusters) <- names(metric_list)
    return(metric_clusters)
}

make_clusters_from_metric_list_mc <- function(metric_list,metric_df=NULL,cutoff_quantile=0.95,metadata_file,perp=10,method="dbscan",n_cores=4,skip_clustering=FALSE){
    require(foreach)
    require(doMC)
    registerDoMC(cores=n_cores)
    source("Andre_F_functions.R")
    

    metric_clusters <- list()
    foreach(i=1:length(metric_list),.combine="c",.multicombine=TRUE) %dopar% {
        print(paste0("Working ",names(metric_list)[i]))
        df <- metric_list[[i]]
        #str(df)

        if(skip_clustering==FALSE){
        tsne_out <- run_tsne(df,perplexity=perp)
        print("tSNE done")

        tsne_out <- get_clusters_tsne(tsne_out,method)

        tsne_out <- get_metadata_tsne(tsne_out, metadata_file)
        #str(rownames(tsne_out))

        #print(head(tsne_out))
        tsne_out <- get_drivers_from_clusters(df,tsne_out)
        #str(tsne_out)
        metric_clusters <- c(metric_clusters, list(tsne_out)); print(paste0("Finished metric ",i))} else if(skip_clustering==TRUE){
            tsne_out <- df
            print("Skipping Clustering")
            tsne_out <- get_drivers_from_clusters(metric_df,tsne_out)
            metric_clusters <- c(metric_clusters, list(tsne_out))
            print(paste0("Finished metric ",i))}
        return(metric_clusters)
    }
 #   names(metric_clusters) <- names(metric_list)

}





make_TF_target_mat <- function(TF,indir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/patientEdgeList/",outdir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/"){
    require(igraph)
    source("/home/anf2034/Andre_F_functions.R")

    target_list <- list()
    print(paste0("Working on ", TF))

file_list <- list.files(indir,"_graph.rds")
#file_list <- file_list[grep("_graph",file_list)]
#file_list <- file_list[1:10]

for(i in file_list){graph <- readRDS(paste0(indir,i))
                    targets <- get_regulatory_targets(graph,TF)
                    target_list <- c(target_list,list(targets))
                    print(paste0("Finished ", grep(i,file_list)," of ",length(file_list)))

                }
saveRDS(target_list, paste0(outdir,TF,"_target_list.rds"))


names(target_list) <- gsub("_graph.rds","", file_list)

union_targets <- c();  for(i in target_list){ union_targets <- unique(c(union_targets,i))}
TF_matrix <- matrix(0,length(target_list),length(union_targets))
rownames(TF_matrix) <- names(target_list)
colnames(TF_matrix) <- union_targets
str(TF_matrix)
    for(i in 1:length(target_list)){ targets <- target_list[[i]]; if(length(targets)> 1){ TF_matrix[i, targets] <- 1} else{ TF_matrix[i,] <- 0}}

target_summary <- colSums(TF_matrix)/nrow(TF_matrix)

print(paste0("Variance for ", TF,"= ",round(var(target_summary),digits=2)))

saveRDS(TF_matrix,paste0(outdir,TF,"_target_matrix.rds"))
saveRDS(target_summary,paste0(outdir,TF,"_target_summary.rds"))


}

make_var_estimate_test <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/var_estimate.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    require(igraph)

    source("/home/anf2034/Andre_F_functions.R")
    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),1]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        str(in_mat)
        sub_mat <- in_mat[in_samples,]
        inter_summary <- round(colSums(in_mat)/nrow(in_mat),digits=2)
        intra_summary <- round(colSums(sub_mat)/nrow(sub_mat),digits=2)

        inter_summary <- 100*table(inter_summary)/sum(table(inter_summary))
        intra_summary <- 100*table(intra_summary)/sum(table(intra_summary))

        subtype_df <- rbind(subtype_df,cbind(inter_summary["1"],intra_summary["1"],sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_variability.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_variability.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}

make_var_estimate <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",metadata_column =2, TF_list="/pbtech_mounts/homes024/anf2034/TF_list_Erica_networks.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_Matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")

    TF_list <- readLines(TF_list)
    if(is.character(metadata_file)){
        metadata <- read.table(metadata_file,header=T,sep=',',stringsAsFactors=F)} else if (is.data.frame(metadata_file)){ metadata <- metadata_file}
#    str(metadata)
   print(sample_type)
    in_samples <- metadata[which(metadata[metadata_column] == sample_type),"Sample"]
    if(length(in_samples) <=1){ print("This cohort is too small for reliable analysis")} else{

    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
#        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]
        inter_summary <- round(colSums(in_mat)/nrow(in_mat),digits=2)
        intra_summary <- round(colSums(sub_mat)/nrow(sub_mat),digits=2)


#        inter_summary <- round(sd(colSums(in_mat))/ncol(in_mat),digits=2)
#        intra_summary <- round(sd(colSums(sub_mat))/ncol(sub_mat),digits=2)

        inter_summary <- 100*table(inter_summary>cutoff)/sum(table(inter_summary>cutoff))
#        print(inter_summary)
        intra_summary <- 100*table(intra_summary>cutoff)/sum(table(intra_summary>cutoff))
#        print(intra_summary)

        subtype_df <- rbind(subtype_df,cbind(inter_summary["TRUE"],intra_summary["TRUE"],sample_type,i),stringsAsFactors=F)
#        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter_type","Intra_type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_variability.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_variability.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}}

}





make_dist_estimate <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/TF_list_Erica_networks.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")

    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),"Sample"]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]


        intra_dist <- sd(dist(sub_mat))/mean(dist(sub_mat))
        inter_dist <- sd(dist(in_mat))/mean(dist(in_mat))


        subtype_df <- rbind(subtype_df,cbind(inter_dist,intra_dist,sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}

make_dist_estimate_test <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",TF_list="/pbtech_mounts/homes024/anf2034/var_estimate.txt",cutoff=0.7,matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_matrices/", return_df=FALSE,outdir="/pbtech_mounts/homes024/anf2034/target_variability/"){
    source("/home/anf2034/Andre_F_functions.R")

    TF_list <- readLines(TF_list)
    metadata <- read.table("TCGA_Pan_Can_metadata.csv",header=T,sep=',',stringsAsFactors=F)
    in_samples <- metadata[which(metadata$cohort == sample_type),"Sample"]
    str(in_samples)
    subtype_df <- data.frame(stringsAsFactors=F)
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        in_samples_real <- intersect(in_samples,rownames(in_mat))
        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]

        intra_dist <- sd(dist(sub_mat))/mean(dist(sub_mat))
        inter_dist <- sd(dist(in_mat))/mean(dist(in_mat))


        subtype_df <- rbind(subtype_df,cbind(inter_dist,intra_dist,sample_type,i),stringsAsFactors=F)
        print(paste0("Done ",i))
    }
    colnames(subtype_df) <- c("Inter-type","Intra-type","Disease_Type","TF")
    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.rds"))
    write.table(subtype_df, paste0(outdir,sample_type,"_target_dist_v2.txt"),quote=F, sep='\t')
    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}

}

gr.mid <- function (x) {
    start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
    return(x)
}

make_graph_from_el <- function(infile, outdir,delete=FALSE, prefix="TCGA"){
    require(igraph)
    el <- as.matrix(read.table(infile, header=F, sep='\t', stringsAsFactors=F))[,1:2]
    filename <- paste0(prefix,unlist(strsplit(infile,"TCGA"))[2])
    graph <- graph_from_edgelist(el)
    saveRDS(graph, paste0(outdir,gsub("edgelist.txt","graph.rds",filename)))
    if(delete==TRUE){ system(paste0("rm ",infile))} else{print("Finished")}

}


rewiring_auc <- function(dir){
    file_list <- list.files(dir)
    file_list <- file_list[grep("_matrix.rds", file_list)]
    out_list <- list()
    for(i in file_list){
        in_mat <- readRDS(paste0(dir,i))
        mat_index_order <- rev(names(sort(rowSums(in_mat))))
        target_all_len <- c()
        target_len <- c()
        blacklist <- c()
        for(j in 1:length(mat_index_order)){
            target_vec <- names(which(in_mat[mat_index_order[j],] == 1))
            target_len <- union(target_len,target_vec)
            target_all_len <- c(target_all_len,length(target_len))


        }
        target_all_len <- as.data.frame(target_all_len,stringsAsFactors=F)
        colnames(target_all_len) <- "Union_Targets"

        target_all_len$Num_Samples <- 1:dim(target_all_len)[1]

        out_list <- c(out_list,list(target_all_len))
#        print(paste0("Finished working ", unlist(strsplit(i,"_"))[1]))
    }
    names(out_list) <- gsub("_target_matrix.rds","",file_list)
    return(out_list)}

rewiring_auc2 <- function(dir){
    file_list <- list.files(dir,"_matrix.rds")
    file_list <- file_list

    out_list <- list()
        for(i in file_list){
        in_mat <- readRDS(paste0(dir,i))
        initial <- names(sort(rowSums(in_mat)))[1]
        blacklist <- c(initial)
        target_len <- names(which(in_mat[initial,] == 1))
        target_all_len <- length(target_len)

        for(j in 2:nrow(in_mat)){
            sub_mat <- in_mat[setdiff(rownames(in_mat),blacklist),]
            if(is.null(nrow(sub_mat))){target_index <- setdiff(rownames(in_mat),blacklist) } else{ sub_mat <- sub_mat
            target_list <- apply(sub_mat, 1, function(x) colnames(sub_mat)[which(x == 1)])
            target_list <- sort(unlist(lapply(target_list, function(x) length(setdiff(x, target_len)))))



            target_index <- names(target_list)[1]}

            blacklist <- c(blacklist,target_index)
            target_vec <- names(which(in_mat[target_index,] == 1))


            target_len <- union(target_len,target_vec)
            target_all_len <- c(target_all_len,length(target_len))



        }
        target_all_len <- as.data.frame(target_all_len,stringsAsFactors=F)
        colnames(target_all_len) <- "Union_Targets"

        target_all_len$Num_Samples <- 1:dim(target_all_len)[1]

        out_list <- c(out_list,list(target_all_len))

        print(paste0("Finished ",gsub("_target_matrix.rds","",i)))





    }
    names(out_list) <- gsub("_target_matrix.rds","",i)
    return(out_list)
    }

rewiring_auc_sequence <- function(dir){
    file_list <- list.files(dir)
    file_list <- file_list[grep("_matrix.rds", file_list)]
#    file_list <- file_list[1:5]
    out_list <- data.frame()
    for(i in file_list){
        in_mat <- readRDS(paste0(dir,i))
        mat_index_order <- rev(names(sort(rowSums(in_mat))))
        target_list <- list()

        for(j in 1:length(mat_index_order)){
            target_vec <- names(which(in_mat[mat_index_order[j],] == 1))
            target_list <- c(target_list,list(target_vec))
        }
        out_list <- rbind(out_list,cbind(gsub("_target_matrix.rds","",i), length(unique(lapply(target_list, sort)))),stringsAsFactors=F)

    }
    out_list[,2] <- as.numeric(out_list[,2])

    return(out_list)}




add_graphs <- function(graph1,graph2){
    if((is.weighted(graph1) & is.weighted(graph2)) == FALSE){
        print("One of the graphs is not weighted, please recheck")} else{
            graph_out <- graph1 + graph2
            E(graph_out)$weight_1[is.na(E(graph_out)$weight_1)] <- 0
            E(graph_out)$weight_2[is.na(E(graph_out)$weight_2)] <- 0
            E(graph_out)$weight <- E(graph_out)$weight_1+E(graph_out)$weight_2
            return(graph_out)}}

subtract_graphs <- function(graph1, graph2){
    if((is.weighted(graph1) & is.weighted(graph2)) == FALSE){
        print("One of the graphs is not weighted, please recheck")} else{
            graph_out <- union(graph1,graph2)
            E(graph_out)$weight_1[is.na(E(graph_out)$weight_1)] <- 0
            E(graph_out)$weight_2[is.na(E(graph_out)$weight_2)] <- 0
            E(graph_out)$weight <- E(graph_out)$weight_1 - E(graph_out)$weight_2
            return(graph_out)}}

read_peak_gene_el <- function(edgelist,cutoff=0.5){

    edgelist <- read.table(edgelist, sep=',',stringsAsFactors=F,header=F)
    split1 <- unlist(strsplit(edgelist[,1],"-(?=[^-]+$)", perl=TRUE))
    edgelist$Peak <- split1[seq(2,length(split1),2)]
    edgelist$Gene <- split1[seq(1,length(split1),2)]
    edgelist <- edgelist[which(edgelist[,2] >cutoff),]
    edgelist <- edgelist[,c(4,5,2)]
    colnames(edgelist) <- c("Peak","Gene","Weight")


    return(edgelist)}

parse_ypredict_outputs <- function(dir, Peak_gr, cutoff=0.5,peak_cutoff=-0.5,peak_height_mat,split="_") {
    require(GenomicRanges)
    file_list <- list.files(dir,".csv")
    for(i in file_list){
        el <- read_peak_gene_el(paste0(dir,i),cutoff)
        if(!is.null(split)){
        sample_ID <- gsub(".csv","",unlist(strsplit(i,"_"))[2])} else{sample_ID <- gsub("yPredict_","",gsub(".csv","",i))}
        
        sample_ID <- gsub("-","", sample_ID)
        str(sample_ID)
        el_int <- Peak_gr[el$Peak]
        gr_names <- el_int$name
        el_int@elementMetadata@listData[setdiff(y="name", names(el_int@elementMetadata@listData))] <- NULL
        el_int$Gene <- el$Gene
        el_int$Weight <- el$Weight
        el_int$score <- peak_height_mat[el_int$name,sample_ID]

        if(is.null(peak_cutoff)){ el_int <- el_int} else{
            el_int <- el_int[which(el_int$score >=peak_cutoff)]
            el_int$score <- NULL
         
        }
        el_int$Weight <- NULL
        gr_to_bed(el_int, gsub(".csv",".bed",paste0(dir,gsub("yPredict_","",i))),TRUE)
        
        print(paste0("Finished ",i))
    }}


weighted_graph_from_edgelist <- function(edgelist, column=3){
    graph_out <- graph_from_edgelist(as.matrix(edgelist[c(1,2)]))
    E(graph_out)$weight <- 1
    return(graph_out)

}

normalize_peak_gene_matrix <- function(matrix,margin='row'){
    if(margin == "row"){
        sample_count <- apply(matrix,1, max)
        matrix_out <- matrix/sample_count
    }
    else if(margin=="column"){
        sample_count <- apply(matrix,2, max)
        matrix_out <- matrix/sample_count
    }
    return(matrix_out)}



#object_size <- function(objects){
#    if(is.vector(objects) == TRUE){



#psetdiff_with_metadata <- function(gr1,gr2){
 #   hits <- findOverlaps}


expand.matrix <- function(A, sparse=TRUE){
    require(Matrix)
    if(sparse ==FALSE){
        m <- nrow(A)
        n <- ncol(A)
        B <- matrix(0,nrow = m, ncol = m)
        C <- matrix(0,nrow = n, ncol = n)
        out <- cbind(rbind(B,t(A)),rbind(A,C))} else{
            m <- nrow(A)
            n <- ncol(A)
            B <- Matrix(0,nrow = m, ncol = m,sparse=TRUE)
            C <- Matrix(0,nrow = n, ncol = n,sparse=TRUE)
            out <- cbind(rbind(B,t(A)),rbind(A,C))}
    return(out)
}

expand.graph.matrix <- function(A,sparse=TRUE){
    require(Matrix)
    if(sparse ==TRUE){
        print("Working 1")
        m <- nrow(A)
        n <- ncol(A)
        B <- Matrix(0,nrow = m, ncol = m,sparse=TRUE)
        C <- Matrix(0,nrow = n, ncol = n,sparse=TRUE)
        int1 <- rbind(B,t(A))
        int2 <- rbind(A,C)
        out <- cbind(int1,int2)} else{
            print("Working 2")
             m <- nrow(A)
        n <- ncol(A)
        B <- matrix(0,nrow = m, ncol = m)
             C <- matrix(0,nrow = n, ncol = n)


             out <- cbind(rbind(B,t(A)),rbind(A,C))}


    out<- triu(out)
    return(out)

}

get_non_zero <- function(mat,margin="row"){
    out_vec <- c()
    for(i in 1:nrow(mat)){
        out_vec <- c(out_vec,all(mat[i,] == 0))
    }
    return(out_vec)

}

get_most_variable <- function(mat, margin=c('row','column'),val_cutoff=NULL,quantile=NULL,verbose=F){
    if(any(is.na(mat)) == TRUE){
        mat[is.na(mat)] <- 0
        if(verbose){print("There are NA values in the matrix. Working anyway")}}

    else{ if(verbose){print("Working")}}


    if(is.null(quantile) & !is.null(val_cutoff)){
        cutoff= val_cutoff
    } else if(is.null(quantile) & is.null(val_cutoff)) { cutoff = 0} else if(!is.null(quantile) & is.null(val_cutoff)){ cutoff = quantile
                                                                   } else { print("Can't have an absolute val_cutoff *AND* quantile-based cutoff; setting to 0"); cutoff =0}

    if(verbose){print(cutoff)}
    if( margin =="row"){
        vec <- apply(mat, 1, var)
        vec <- sort(vec,decreasing=TRUE)
  
    }
    else if( margin== "column"){
        vec <- apply(mat,2,var)
    }

    if(!is.null(val_cutoff)){ variable <- which(vec >= cutoff)} else { variable <- which(vec > quantile(vec,cutoff,na.rm=TRUE))}

    return(variable)}

sparse_mat_to_vector <- function(mat){
    mat <- as(mat, "dgTMatrix")
    out <- mat@x
    names(out) <- paste0(rownames(mat)[(mat@i+1)],"~",colnames(mat)[(mat@j+1)])

    return(out)}

cutdown.graph.matrix <- function(graph_mat){
    graph_sub_mat <- graph_mat[grep("[A-Z]_[0-9]", rownames(graph_mat)),]
    graph_sub_mat <- graph_sub_mat[,setdiff(colnames(graph_sub_mat),rownames(graph_sub_mat))]
    return(graph_sub_mat)}


make_combined_metric <- function(outdegree_df,expression_df, betweenness_df=NULL,method=c("simple","complex")){
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


get_link_abundance <- function(P_G_link_df, metadata_df, column){
    require(dplyr)
    sample_groups <- lapply(unique(metadata_df[,column]), function(x) dplyr::filter(metadata_df, get(column)==x)$Sample)
    names(sample_groups) <- unique(metadata_df[,column])
    sample_groups <- sample_groups[unlist(lapply(sample_groups, length) >1)]

    abundance_mat <- do.call("cbind",lapply(sample_groups, function(x) rowSums(P_G_link_df[,x])/length(x)))
    colnames(abundance_mat) <- names(sample_groups)
    return(abundance_mat)
}


get_tissue_specific_links <- function(abundance_mat,zscore=2.5, min_fraction=0.5, output_format=c("original","zscore","both")){
    sub <- zscore_df(abundance_mat, "row")
    zscore_mat <- sub
    sub <- apply(sub,1, function(x) any(x >=zscore))
    sub <- names(which(sub == TRUE))
    print(paste0(length(sub)," links survived first cut"))
    out_mat <- abundance_mat[sub,]
    index <- apply(out_mat,1, function(x) any(x >=min_fraction))
    index <- index[which(index == TRUE)]
    print(paste0(length(index)," links survived second cut"))
    out_mat <- out_mat[names(index),]
    out_mat <- as.data.frame(out_mat)

    zscore_mat <- as.data.frame(zscore_mat[names(index),])
    if(output_format=="original"){ final <- out_mat} else if(output_format=="zscore"){ final <- zscore_mat} else if(output_format=="both") { final <- list(out_mat,zscore_mat); names(final) <- c("Original","Zscore")}
        
    return(final)
}


get_coordinate <- function(peak,gr,peak_column){
    metadata <- unlist(gr@elementMetadata@listData[peak_column])
    index <- grep(peak,metadata)
    chrom <- as.character(gr@seqnames)[index]
    start <- gr@ranges@start[index]
    end <- start+gr@ranges@width[index]
    Peaks <- peak
    coord_df <- cbind(chrom,start,end,Peaks)
    return(coord_df)}

PWM_to_Jaspar_norm <- function(PWM){
    PWM_name <- names(PWM@name)

    PWM_mat <- PWM@profileMatrix
    PWM_out <- list(paste0(">",PWM_name))
    for(i in 1:nrow(PWM_mat)){ vec <- paste0(rownames(PWM_mat)[i]," \t"," [ ", paste0(unlist(PWM_mat[i,]),collapse=" ")," ] ")
                               PWM_out <- c(PWM_out,list(vec))}
    return(PWM_out)
}

PIQ_bed_to_gr <- function(bed){
require(GenomicRanges)
    check <- readLines(bed)
    search_check <- grep("track", check)
    if(length(search_check) >0){
       new_bed <- gsub(".bed", "_fixed.bed",bed)
    system(paste0("sed -i '/^track/d' ",bed))

        #bed <- new_bed
        print("Bed file was incorrectly formatted")} else{ print("Bed file was correctly formatted")}
     bed <- read.table(bed, stringsAsFactors=F, sep="\t",header=F)
                                 gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),purity=bed$V5,strand=bed$V6,motif=bed$V4)
    return(gr)}

PIQ_bed_to_gr_v2 <- function(bed){
    require(GenomicRanges)
    system(paste0("sed -i '/^track/d' ",bed))
    bed <- read.table(bed, stringsAsFactors=F, sep="\t",header=F)
    gr <- GRanges(seqnames=Rle(bed$V1), ranges= IRanges(start=bed$V2, end=bed$V3),purity=bed$V5,strand=bed$V6,motif=bed$V4)
    return(gr)}


get_real_peaks <- function(peak_patient_mat,cutoff=0.05,output="names"){
    peak_list <- list()

    for(i in 1:ncol(peak_patient_mat)){
        sub <- as.data.frame(peak_patient_mat[,i],stringsAsFactors=F)
        colnames(sub) <- "sub"
        rownames(sub) <- rownames(peak_patient_mat)
        sub$quantile <- ecdf(sub$sub)(sub$sub)
        sub$pval <- 1-sub$quantile
        sub$Peak <- rownames(sub)
        real_peaks <- sub[which(sub$pval <= cutoff),]

        if(output == "names"){
            peaks_out <- rownames(real_peaks)
        } else if (output != "names"){
            peaks_out <- real_peaks}

        peak_list[[colnames(peak_patient_mat)[i]]] <- peaks_out}
    return(peak_list)
}


get_real_peaks_v2 <- function(peak_patient_mat,cutoff=1e-5,output="names"){
     peak_list <- list()

    for(i in 1:ncol(peak_patient_mat)){
        sub <- as.data.frame(peak_patient_mat[,i],stringsAsFactors=F)
        colnames(sub) <- "sub"
        rownames(sub) <- rownames(peak_patient_mat)
        sub$pval <- 1-ppois(2^sub$sub,median(2^sub$sub))
        real_peaks <- sub[which(sub$pval <= cutoff),]

        if(output == "names"){            peaks_out <- rownames(real_peaks)
        } else {
            peaks_out <- real_peaks}

        peak_list[[colnames(peak_patient_mat)[i]]] <- peaks_out}
    return(peak_list)
}

rank_mat <- function(matrix,margin="column",convert_NA=TRUE){
    if(any(is.na(matrix)) == TRUE){
        if(convert_NA==TRUE){
        matrix[is.na(matrix)] <- 0
        print("There are NA values in your matrix, converting to 0")} else {print("Not converting NA values, will rank them last")}} else{ matrix= matrix}

    if(margin == "column"){
        mat <- apply(matrix,2, function(x) (length(x)+1)-rank(x,na.last=TRUE))}
    else {         mat <- t(apply(matrix,1, function(x) (length(x)+1)-rank(x,na.last=TRUE)))}
    mat <- as.data.frame(mat)
    mat[is.na(mat)] <- 0
    return(mat)}

get_unfinished_PIQ <- function(PIQ_dir,required_files){
    unfinished_list <- list()
    file_list <- list.files(PIQ_dir)
    for(i in file_list) {
        indir <- paste0(PIQ_dir,i)
        present <- list.files(indir, ".bed")
        missing <- setdiff(required_files,present)
        unfinished_list[[i]] <- missing}

    return(unfinished_list)}

parse_missing <- function(unfinished_list,motif_file_dir){
    require(stringr)
    out_list <- list()
    for(i in 1:length(unfinished_list)){
        sub <- unfinished_list[[i]]
        if(length(sub) >0){

        sample <- names(unfinished_list)[i]
        sub2 <- lapply(sub, function(x) unlist(strsplit(x,"-")))
        sub_names <- lapply(sub2, function(x) x[2])
        sub_names <- unique(gsub(".RC","",sub_names))


        search_string1 <- lapply(sub_names, function(x) str_extract(x,"LINE[0-9]+"))
        search_string2 <- lapply(sub_names, function(x) unlist(strsplit(x,"LINE[0-9]+"))[1])



        motif_file1 <- lapply(search_string1, function(x) system(paste0("grep -r -l ", x," ",motif_file_dir),intern=T))
        motif_file2 <- lapply(search_string2, function(x) system(paste0("grep -r -l ", x," ",motif_file_dir),intern=T))

#        motif_file <- lapply(1:length(motif_file1), function(x) intersect(motif_file1[[x]],motif_file2[[x]]))
        motif_file <- lapply(1:length(motif_file1), function(x) c(motif_file1[[x]],motif_file2[[x]]))



 #       str(motif_file1)
  #      str(motif_file2)
  #      str(search_string1)
  #      str(search_string2)
  #      str(motif_file)
        motif_file <- lapply(motif_file, function(x) names(which(table(x) == max(table(x)))))




        out_df <- as.data.frame(cbind(sample,unlist(search_string1),unlist(motif_file)),stringsAsFactors=F)
        colnames(out_df) <- c("Sample","Motif_Symbol","Motif_File")
        out_list[[sample]] <- out_df} else{ print("No missing TFs")}



    }
        return(out_list)
}

make_PIQ_step1_job_table <- function(PIQ_table, parsed_missing_list,ext="_ATAC_hg38"){
    out_table <- data.table()

    for(i in 1:length(parsed_missing_list)){
        sub <- parsed_missing_list[[i]]
        DHS_file <- grep(names(parsed_missing_list)[i], PIQ_table$DHS_file,value=T)
        Motif_Symbol <- sub[,2]
        Motif_file <- sub[,3]
        tmp_folder <- gsub(".txt","",unlist(lapply(sub[,3], function(x) unlist(strsplit(x,"_PWM_"))[2])))
        ID <- gsub(ext,"",names(parsed_missing_list)[i])
        ID <- paste0(ID,"_cleanup_",1:length(tmp_folder))

        sub_table <- data.table(DHS_file=unique(DHS_file),Motif_File=Motif_file,PIQ_Output_Dir=unique(PIQ_table$PIQ_Output_Dir),Final_Output_Dir=unique(PIQ_table$Final_Output_Dir), tmp_folder=tmp_folder,ID=ID,Motif_Symbol=Motif_Symbol)
        out_table <- rbind(out_table,sub_table)}
    return(out_table)}

CCLE_Copy_Number_to_gr <- function(CCLE_csv_file){
    require(GenomicRanges)
    require(stringr)
    infile <- read.table(CCLE_csv_file,header=T,sep=",", stringsAsFactors=F)
    infile$Type <- unlist(lapply(infile$sample, function(x) unlist(str_split(x,"_",n=2))[2]))
    infile$Cell_Line <- unlist(lapply(infile$sample, function(x) unlist(str_split(x,"_",n=2))[1]))

    gr <- GRanges(seqnames=Rle(infile[,"Chromosome"]), ranges= IRanges(start=infile[,"Start"], end=infile[,"End"]),strand="*",Cell_Type=infile[,"Type"],Cell_Line=infile[,"Cell_Line"],Modal_CN=infile[,"Modal_Total_CN"],Homozyg_Del=infile[,"Homozygous_deletion"],Num_Probes=infile[,"Num_Probes"],Region_Length=infile[,"Length"])
    return(gr)}

plot_peak_gene_link_CN <- function(CCLE_gr,peak_gene_list,peak_gr){
    require(GenomicRanges)
    require(ggplot2)
    peaks_in <- unlist(lapply(peak_gene_list, function(x) unlist(strsplit(x,"~"))[1]))
    str(peaks_in)
    int_gr <- intersect_with_metadata(CCLE_gr,peak_gr[peaks_in])
    found_peaks <- int_gr[unique(findOverlaps(CCLE_gr,int_gr)@to)]
    str(found_peaks)
    if(length(int_gr) >0){
        int_df <- data.frame(int_gr$Cell_Line,int_gr$Modal_CN,int_gr$Num_Probes,peak_gene_list,found_peaks)
        colnames(int_df)[4:5] <- c("Links","Peaks")
        colnames(int_df) <- gsub("int_gr.","",colnames(int_df))
        str(int_df)
        out <- ggplot(int_df,aes(Cell_Line,Links, fill=Modal_CN))+geom_raster()+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=8),axis.text.y=element_text(face='bold',size=15))
                                        #    print(int_df)
    } else{ print("No copy number information for these regions")}

    return(out)}

plot_region_CN <- function(CCLE_gr,peak_gene_list,peak_gr){
    require(GenomicRanges)
    require(ggplot2)
    peaks_in <- unlist(lapply(peak_gene_list, function(x) unlist(strsplit(x,"~"))[1]))
    str(peaks_in)
    int_gr <- intersect_with_metadata(CCLE_gr,peak_gr[peaks_in])
    if(length(int_gr) >0){
        int_df <- data.frame(int_gr$Cell_Line,int_gr$Modal_CN,int_gr$Num_Probes,peak_gene_list,peaks_in)
        colnames(int_df)[4:5] <- c("Links","Peaks")
        colnames(int_df) <- gsub("int_gr.","",colnames(int_df))
        out <- ggplot(int_df,aes(Cell_Line,Links, fill=Modal_CN))+geom_raster()+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=8),axis.text.y=element_text(face='bold',size=15))
                                        #    print(int_df)
    } else{ print("No copy number information for these regions")}

    return(out)}



uncompress_PIQ <- function(PIQ_gz_path,outdir){
    system(paste0("mkdir ", outdir))
    outdir <- paste0(outdir,"/")

    infile <- PIQ_gz_path
    outfile <- gsub(".gz","",PIQ_gz_path)

    if(file.exists(outfile) == FALSE){

        system(paste0("gunzip -c ",infile," > ",outfile))
    } else if(file.exists(outfile)==TRUE){ print("Not decompressing file")}

    system(paste0("csplit -f ",outdir,"/out ", outfile," '/track/' '{*}' --suppress-matched -z"))
    file_list <- list.files(outdir,"out[0-9]")
    str(file_list)
    for(i in file_list){
        int_file <- paste0(outdir,i)
        str(int_file)
        name <- system(paste0("awk '{print $4 }' ",outdir,i,"|head|uniq"),intern=T)
        str(name)
        system(paste0("mv ",outdir,i," ",outdir,name,"-calls.all.bed"))}

    system(paste0("rm ",outfile))

}

calc_metric <- function(indir, outfile,metric=c("outdegree","indegree","betweenness")){
    require(igraph)

    file_list <- list.files(indir,"_graph.rds")

    metric_df <- data.frame(0)
    rownames(metric_df) <- "Empty"
    str(file_list)

    print(paste0("Calculating ",metric))

    for(i in file_list){
#        str(i)
        graph <- readRDS(paste0(indir,i))

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

    print(paste0("Finished processing ", grep(i,file_list)," of ", length(file_list)))
    }
    colnames(metric_df) <- gsub("_Final|_PIQ|_graph.rds","",colnames(metric_df))
    metric_df <- metric_df[-1,-1]
    metric_df[is.na(metric_df)] <- 0

saveRDS(metric_df,outfile)
}


outdegree <- function(graph){
    require(igraph)
    out <- degree(graph,mode='out')
    out <- out[which(out !=0)]
    out <- out[sort(names(out))]
    return(out)}

subtract_dfs <- function(df1,df2){
    common_rows <- intersect(rownames(df1),rownames(df2))
    common_cols <- intersect(colnames(df1),colnames(df2))

    diff_df <- df1[common_rows,common_cols]-df2[common_rows,common_cols]
    return(diff_df)}

get_non_redundant_columns <- function(mat,cutoff=0.99,TF_margin=c("row","column")){
    require(Hmisc)
    mat_back <- mat

    if(TF_margin == "row"){
        mat <- rcorr(t(as.matrix(mat)))}
    else{ mat <- rcorr(as.matrix(mat))}


    mat <- mat[[1]]
##    str(mat)
    mat[is.na(mat)] <- 1

    mat[upper.tri(mat)] <- 0
    diag(mat) <- 0
##    print(mat[1:10,1:10])
    mat_out <- mat[,!apply(mat,2,function(x) any(x >= cutoff))]

    mat_out <- mat_back[colnames(mat_out),]
##    str(mat_out)
    return(mat_out)
}

compare_graphs <- function(indir){
    require(igraph)
    require(Matrix)

    file_list <- list.files(indir,"_graph.rds")
#    file_list <- file_list[1:5]
    out_list <- list()
    print("Making list of edges by sample ")

    for(i in file_list){
        sample_name <- gsub("_Final|_PIQ|_graph.rds","",i)
        graph <- readRDS(paste0(indir,i))
        graph_mat <- get.adjacency(graph)
        vec <- sparse_mat_to_vector(graph_mat)
        out_list[sample_name] <- list(vec)
    }
    print("Finished making list of edges for all samples")

    edge_names <- unique(unlist(lapply(out_list,names)))
    mat <- Matrix(0,nrow=length(edge_names),ncol=length(names(out_list)), dimnames=list(edge_names,names(out_list)),sparse=TRUE)
    str(mat)
    str(edge_names)
    print("Making giant matrix, hold on to your pants")
    for(i in names(out_list)){print(paste0("Working ",i)); sub <- out_list[[i]]; index <- names(sub); mat[index,i] <- 1; print(paste0("Finished ",grep(i,names(out_list))," of ", length(names(out_list))))}

    final_out <- list(out_list,mat)

    return(final_out)}

compare_graphs2 <- function(indir){
    require(igraph)
    require(Matrix)

    file_list <- list.files(indir,"_graph.rds")
    #file_list <- file_list[1:5]
    out_list <- list()
    print("Making list of edges by sample ")

    for(i in file_list){
        sample_name <- gsub("_Final|_PIQ|_graph.rds","",i)
        graph <- readRDS(paste0(indir,i))
        graph_mat <- get.adjacency(graph)
        vec <- sparse_mat_to_vector(graph_mat)
        out_list[sample_name] <- list(vec)
    }
    print("Finished making list of edges for all samples")

    edge_names <- lapply(out_list,function(x) names(x))
    edge_names <- do.call(rbind,lapply(edge_names,function(x) as.data.frame(cbind(gsub("_Final|_graph.rds|_PIQ","",names(edge_names)),x),stringsAsFactors=F)))
    edge_names$vals <- 1



    print("Making giant matrix, hold on to your pants")



    rows <- as.numeric(as.factor(edge_names[,2]))
    cols <- as.numeric(as.factor(edge_names[,1]))
    vals <- edge_names[,3]



#    mat <- Matrix(0,length(unique(rows)),length(unique(cols)),dimnames=list(rows,cols))

    mat <- sparseMatrix(i=rows,j=cols,x=vals,dimnames=list(unique(edge_names[,2]),unique(edge_names[,1])))



    final_out <- list(mat,out_list,edge_names)

    return(final_out)}




matrix_density <- function(matrix){

    if(is(matrix,"sparseMatrix") == TRUE){
        dens_out <- length(matrix@x)/(dim(matrix)[1]*dim(matrix)[2])
    } else { print("This will only work for sparse matrices")}
    return(dens_out)
}

analyzing_giant_graph_matrix <- function(matrix, chunk_size=10000,margin="column"){
    #matrix <- matrix[which(rowSums(matrix) <341),]
    blocks <- (nrow(matrix)%/% chunk_size)
    if((nrow(matrix)%% chunk_size) == 0){
        blocks <- blocks} else {blocks <- blocks+1}
    str(blocks)
    #str(matrix)
    if(margin == "column"){
        dist_mat <- matrix(0,nrow=ncol(matrix),ncol=ncol(matrix))
        for(i in 1:blocks){
            start <- ((i-1)*chunk_size)+1
            end <- i*chunk_size
            if(end > nrow(matrix)){ end <- nrow(matrix)}
            else{ end <- end}
            sub_mat <- matrix[start:end,]
            nrow(sub_mat)
            if(is.null(nrow(sub_mat)) == FALSE){
                sub_mat <- t(matrix[start:end,])} else {sub_mat <- sub_mat}

            out <- as.matrix(dist(sub_mat))

#            str(out)
#            dist_mat <- dist_mat+out
            dist_mat <- dist_mat+zscore_mat(out)
            #print(summary(as.vector(dist_mat)))

            print(paste0("Working ",i," block of ",blocks))
        }


    }
    else if(margin == "row"){
        dist_mat <- matrix(0,nrow=nrow(matrix),ncol=nrow(matrix))
        for(i in 1:blocks){
            start <- ((i-1)*chunk_size)+1
            end <- i*chunk_size
            sub_mat <- matrix[start:end,]
            out <- as.matrix(dist(sub_mat))
#           dist_mat <- dist_mat+out
            dist_mat <- dist_mat+zscore_mat(out)
            print(paste0("Working ",i," block of ",blocks))
        }



    }
    return(dist_mat)}


zscore_mat <- function(matrix){

    out <- (matrix-mean(as.vector(matrix)))/sd(as.vector(matrix))

    return(out)}

normalizing_giant_mat <-  function(giant_mat,expression_df){
    mat_edge_TFs <- unique(unlist(lapply(rownames(giant_mat), function(x) unlist(strsplit(x,"~"))[1])))
    TF_list <- intersect(mat_edge_TFs,rownames(expression_df))
    sample_list <- intersect(colnames(giant_mat), colnames(expression_df))
    expression_sub <- expression_df[TF_list,]
    expression_sub <- expression_sub[,sample_list]




    for(i in TF_list){
        index <- grep(i,rownames(giant_mat))
        sub_mat <- (giant_mat[index,sample_list]*(unlist(1-log10(expression_sub[i,sample_list]))))
        giant_mat[index,sample_list]@x <- sub_mat@x

        print(paste0("Finished ", i))

    }
    return(giant_mat)
}



normalizing_giant_mat2 <- function(giant_mat, expression_df){
    print("Getting TFs from Matrix")

    mat_edge_TFs <- unique(unlist(lapply(rownames(giant_mat), function(x) unlist(strsplit(x,"~"))[1])))

    TF_list <- sort(intersect(mat_edge_TFs,rownames(expression_df)))

    print("Creating compatible objects")
    sample_list <- intersect(colnames(giant_mat), colnames(expression_df))
    expression_sub <- expression_df[TF_list,]
    expression_sub <- expression_sub[,sample_list]

    print("Getting indices")
    row_indices <- lapply(TF_list, function(x) grep(paste0(x,"~"), rownames(giant_mat)))
#    str(row_indices)
    print("Making sub matrices")

    sub_mat_list <- lapply(row_indices, function(x) giant_mat[x,sample_list])
#    str(sub_mat_list)


    expression_mat_list <- lapply(TF_list, function(x) unlist(expression_sub[x,]))

#    str(expression_mat_list)

    print("Performing matrix multiplication")
    norm_mat_list <- lapply(1:length(sub_mat_list), function(x) t(t(sub_mat_list[[x]]) *round((1/expression_mat_list[[x]]),digits=4)))
    saveRDS(norm_mat_list,"norm_mat_list.rds")

    print("Finishing up")
    len_vector <- 1:length(norm_mat_list)
    split_vector <- split_vector(len_vector,300)

    int_mat <- lapply(split_vector, function(x) do.call("rbind",norm_mat_list[x]))
    out_mat <- do.call("rbind", int_mat)

#    out_mat <- Reduce(rbind,norm_mat_list[-1], norm_mat_list[[1]])
    return(out_mat)

}

get_Achilles_disease_submatrix <- function(matrix,metadata_file ,term, margin=c("row","column")){
    index <- grep(term,metadata_file[,2])
    if(length(index) >0){
        if(margin == "row"){
            out <- matrix[index,]
        }
        else if (margin == "column"){

            out <- matrix[,index]
        }
    } else { print(paste0("Can't find ", term))}



    return(out)}


gg_achilles_boxplot <- function(df,metadata_df,subset=NULL,highlight=NULL,threshold=-0.5) {
    require(ggplot2)
    require(reshape2)



    metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(df)),]
    df <- df[,intersect(colnames(df),metadata_df[,1])]


    if(is.null(subset)) {
        print("Working with all the data")



    } else{

            index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
            df <- df[,index]
        }


    df <- as.data.frame(t(df))
    df$Disease <- metadata_df[rownames(df),grep("TCGA",colnames(metadata_df))]
    df$Cell_Line <- rownames(df)

    melted_mat <- reshape2::melt(df)
    melted_mat$Essential <- melted_mat$value <= threshold
    saveRDS(melted_mat,"test_df.rds")

#    for(i in 1:2){ melted_mat[,i] <- as.character(melted_mat[,i])}

#    str(melted_mat)

    if(is.null(highlight)){
        plot <- ggplot(melted_mat, aes(Disease, y=value, fill=Disease))+geom_violin(aes(alpha=0.6))+geom_boxplot(width=0.05,fill="white")+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=10),strip.text.y = element_text(colour = "black", face = "bold",size=10))+facet_grid(variable~.)+guides(alpha=FALSE)+geom_hline(yintercept=0,linetype=2, color="red",alpha=0.5)+ylim(-1,1)} else { print(paste0("Trying to highlight ",highlight))

                                                                                                                                                                                                                                                                                                       index <- unlist(lapply(highlight, function(x) grep(x,melted_mat$Disease)))
                                                                                                                                                                                                                                                                                                        true_index <- setdiff(1:nrow(melted_mat),index)
                                                                                                                                                                                                                                                                                                        melted_mat$Disease[true_index] <- "Other"
                                                                                                                                                                                                                                                                                                                plot <- ggplot(melted_mat, aes(Disease, y=value, fill=Disease))+geom_violin(aes(alpha=0.6))+geom_boxplot(width=0.05,fill="white")+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=10),strip.text.y = element_text(colour = "black", face = "bold",size=10))+facet_grid(variable~.)+guides(alpha=FALSE)+geom_hline(yintercept=0,linetype=2, color="red",alpha=0.5)+ylim(-1,1)}

    return(plot)}

achilles_heatmap <- function(df,metadata_df,subset=NULL,color_palette_list,all_colors=FALSE,title=NULL,scale=c("none","row","column")) {
    require(gplots)
    require(RColorBrewer)

    metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(df)),]
    df <- df[,intersect(colnames(df),metadata_df[,1])]


    if(is.null(subset)) {
        print("Working with all the data")
        color_list <- color_palette_list




    } else{
        if (all_colors == FALSE){
            color_list <- color_palette_list[intersect(names(color_palette_list),subset)]
            sub_col <- color_palette_list[setdiff(names(color_palette_list),subset)]
            sub_col <- lapply(sub_col, function(x) x="gray74")
            color_list <- c(color_list, sub_col)} else {color_list <- color_palette_list}

        #index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
        #df <- df[,index]
    }
    mat2 <- as.matrix(df)

    colors <- unlist(color_list[unlist(metadata_df[colnames(mat2),"TCGA"])])

    p <- try(heatmap.2(mat2,trace='none',cexRow=1.5, cexCol=0.9,margins=c(6,8),symkey=FALSE, symbreaks=TRUE, col=brewer.pal(7,"RdBu"),dendrogram="none",density.info="none",keysize=1,key.title=FALSE,key.xlab="Proliferation rate after \n knockdown",scale=scale,breaks=seq(-1,1,length.out=8), ColSideColors=colors,main=title))
    legend(x="bottomleft",y,title="TCGA equivalent",legend=names(color_list),fill=unlist(color_list),cex=1.2)



#  return(p)
}

achilles_heatmap_binary <- function(df,metadata_df,subset=NULL,color_palette_list,all_colors=FALSE,title=NULL,threshold=-0.75) {
    require(gplots)
        require(gplots)

    metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(df)),]
    df <- df[,intersect(colnames(df),metadata_df[,1])]
    df <- df <= threshold

    df2 <- apply(df,2, function(x) as.numeric(x))
    rownames(df2) <- rownames(df)
    df <- df2
#    str(df)


    if(is.null(subset)) {
        print("Working with all the data")
        color_list <- color_palette_list




    } else{
        if (all_colors == FALSE){
            color_list <- color_palette_list[intersect(names(color_palette_list),subset)]
            sub_col <- color_palette_list[setdiff(names(color_palette_list),subset)]
            sub_col <- lapply(sub_col, function(x) x="gray74")
            color_list <- c(color_list, sub_col)} else {color_list <- color_palette_list}

#        str(color_list)
        #index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
        #df <- df[,index]
    }
    mat2 <- as.matrix(df)

    colors <- unlist(color_list[unlist(metadata_df[colnames(mat2),"TCGA"])])

    p <- try(heatmap.2(mat2,trace='none',cexRow=1.5, cexCol=0.9,margins=c(6,8),symkey=FALSE, symbreaks=TRUE, col=c("white","gray27"),dendrogram="none",density.info="none",keysize=1,key.title=FALSE,key.xlab="Proliferation rate after \n knockdown",scale="none", ColSideColors=colors,main=title))
    legend(x="bottomleft",y,title="TCGA equivalent",legend=names(color_list),fill=unlist(color_list),cex=1.2)




}



achilles_comp <- function(df, metadata_df,driver_list, threshold=0.1){
    out_df <- data.frame(stringsAsFactors=F)
    metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(df)),]
    df <- df[,intersect(colnames(df),metadata_df[,1])]

    for(i in driver_list){
        driver_set <- unlist(strsplit(i[,1],"-"))
        driver_set <- intersect(driver_set, rownames(df))
        disease <- unlist(strsplit(i[,3],"-"))
        cluster <- i[,2]
        gene_set <- rownames(df)

        in_index <- grep(paste0(disease,collapse="|"),metadata_df[,"TCGA"])
        out_index <- setdiff(1:nrow(metadata_df),in_index)

        in_avgs <- unlist(lapply(gene_set, function(x) mean(unlist(df[x,in_index]),na.rm=TRUE)))
        out_avgs <- unlist(lapply(gene_set, function(x) mean(unlist(df[x,out_index]),na.rm=TRUE)))

        in_sd <- unlist(lapply(gene_set, function(x) sd(unlist(df[x,in_index]),na.rm=TRUE)))
        out_sd <- unlist(lapply(gene_set, function(x) sd(unlist(df[x,out_index]),na.rm=TRUE)))

        if(length(in_index ) <2){ in_avgs <- in_avgs
                              } else { in_avgs <- in_avgs*in_sd}

        out_avgs <- out_avgs*out_sd


        if(all(is.na(in_avgs)) == FALSE){
            pvals <- unlist(lapply(gene_set, function(x) wilcox.test(unlist(df[x,in_index]),unlist(df[x,out_index]))$p.value))} else{ pvals <- rep(1,length(out_avgs))}

        disease_out <- paste0(disease,collapse="-")
        sub_df <- as.data.frame(cbind(cluster,in_avgs,out_avgs,pvals,disease_out,gene_set),stringsAsFactors=F)
        sub_df$Driver <- "No"
        sub_df[which(sub_df[,6] %in% driver_set),"Driver"] <- "Yes"
        sub_df$Q_value <- p.adjust(sub_df$pvals,method="BH")
        out_df <- rbind(out_df, sub_df,stringsAsFactors=F)
#        str(sub_df)
        print(paste0("Finished ", i[,3]))


    }
    colnames(out_df) <- c("Cluster","Cluster_Avg","Out_Avg","P_value","TCGA_Disease","Gene","Driver","Q_value")

    for(i in 2:4){ out_df[,i] <- as.numeric(out_df[,i])}

    out_df$Significant <- "No"
    out_df[which(out_df$Q_value <= threshold),"Significant"] <- "Yes"
    out_df$Label <- paste0(out_df$Cluster,":",out_df$TCGA_Disease)
    out_df$Eval <- paste0(out_df$Driver,"-",out_df$Significant)



    return(out_df)}

make_driver_list <- function(metric_cluster_df,column=c("Cluster_drivers","Signif_drivers")){
    metric_cluster_df$Cluster_Majority <- "None"
    for(i in unique(metric_cluster_df$Cluster)){
        index <- which(metric_cluster_df$Cluster == i)
        summary <- table(metric_cluster_df[index,"Subtype"])

        metric_cluster_df_maj <- names(which(summary == max(summary)))
        if(length(metric_cluster_df_maj) > 1){
            metric_cluster_df_maj <- metric_cluster_df_maj[1]
            print("Multiple disease types are tied")}
        else{ }

        metric_cluster_df[index,"Cluster_Majority"] <- metric_cluster_df_maj
    }

    driver_list <- list()
    for(k in unique(metric_cluster_df$Cluster)){
        index <- which(metric_cluster_df$Cluster == k)
        drivers <- unique(metric_cluster_df[index,column])
        subtype_table <- table(metric_cluster_df[index,"Subtype"])
        fraction <- subtype_table/sum(subtype_table)
        fraction <- fraction[which(fraction >=0.05)]
#        str(fraction)


        metric_cluster_subtypes <- paste0(names(fraction),collapse="-")
        clust_id <- paste0("Cluster-",k)
        driver_list[[clust_id]] <- as.data.frame(cbind(drivers,clust_id,metric_cluster_subtypes),stringsAsFactors=F)


    }
    return(driver_list)}



plot_regions_track <- function(plot_ranges, track_names=NULL ,encode_track=NULL, input=NULL,y.field=NULL,buffer=1e5,lines=FALSE,track_smooth=FALSE,smooth_n=50,track_col=NA,bars=FALSE,links=NULL,gencode_collapse=TRUE,hide_names=TRUE,y_cap=NULL,y0=NA,y1=NA,plot_name="",y_height=1,normalize_height=TRUE,...){
    require(gTrack)
    require(gUtils)
    require(data.table)
    require(GenomicRanges)


    source("/home/forbesa1/Andre_F_functions.R")

    if(is.list(input)){ all_ranges <- input} else{
        all_ranges <- list(...)}
    if(length(track_names) != length(all_ranges)){
        print("The number of tracks is not equal to the number of names provided")
    } else { for(i in 1:length(all_ranges)){
                 all_ranges[[i]]$Track_name <- track_names[i]
                 if(hide_names==TRUE){ names(all_ranges[[i]]) <- NULL} else{ names(all_ranges[[i]]) <- names(all_ranges[[i]])}}}
         
    ##print(unlist(lapply(all_ranges, function(x) unique(x$Track_name))))
    ##print(plot_ranges)
    all_ranges <- lapply(all_ranges, function(x) unique(intersect_with_metadata(x,plot_ranges+buffer)))
    names(all_ranges) <- track_names
    if(length(track_col)==1){ track_col <- rep(track_col,length(all_ranges))} else if(length(track_col)==length(all_ranges)){ track_col <- track_col} else{ track_col <- sample(track_col,length(all_ranges))}

    ##print(track_col)

    if(normalize_height==TRUE){
        if(is.list(input) & !is.null(y.field)){
        ymax <- max(unlist(lapply(input, function(x) max(mcols(intersect_with_metadata(x,plot_ranges))[,y.field])))); y1=ymax} else if (!is.null(y.field)){ ymax <- max(mcols(input)[,y.field]); y1 <- ymax}}

    if(is.null(y.field)){
        all_gt <- lapply(1:length(all_ranges), function(x) gTrack(all_ranges[[x]],col=track_col[x],name=unique(all_ranges[[x]]$Track_name,...)))} else{
           all_gt <- list()
            if(track_smooth== TRUE){
                for(x in names(all_ranges)){
                    sub_gr <- all_ranges[[x]]
##                    print(sub_gr)
                    sub_gr <- smooth_gr(sub_gr,binsize=smooth_n,pad=smooth_n)
                    all_ranges[[x]] <- sub_gr
  ##                  roll_mean <- data.table::frollmean(values(all_ranges[[x]])[,y.field],smooth_n,fill=0)
                    #print(roll_mean)
                    #print(values(all_ranges[[x]])[,y.field])

##                    values(all_ranges[[x]])[,y.field] <- roll_mean

                }} else{ print("Not smoothing intervals")}
##            print(all_ranges)
            index <- unlist(lapply(all_ranges,function(x) grep(y.field,names(x@elementMetadata@listData))))
##            str(index)

##            print(length(all_ranges))
##            print(length(track_names))

            for(i in 1:length(index)){
  ##              print(i)
                if(index[i] != 0){
                    ##print("Working 1")
                    sub_gr <- all_ranges[[i]]
  ##                  print(sub_gr)
  ##                  print(names(index)[i])
                    if(length(grep("^Peaks|^Anchors",unique(sub_gr$Track_name))) ==0){
        ##                print("Sub 1")
                        ##  print(track_col[i])
                        
                        sub_gt <- gTrack(sub_gr,name=track_names[i],y.field=y.field,bars=bars,lines=lines,col=track_col[i],y0=y0,y1=y1,...)}  else if(length(grep("^Peaks|^Anchors",track_names[i])) >0) {
      ##                        print("Sub 2")
                            sub_gt <- gTrack(sub_gr,name=track_names[i],y.field=y.field,lines=FALSE,y0=y0,y1=y1,...)}
                                   all_gt[[track_names[i]]] <- sub_gt
                } else{ ##print("Working 2")
                    sub_gr <- all_ranges[[i]]
                                       sub_gt <- gTrack(sub_gr, name=track_names[i],lines=lines,...)
                                                    all_gt[[track_names[i]]] <- sub_gt
                                   }}}

#    print(all_gt)


    print("Collapsing")
    out_track <- collapse_gtrack_list(all_gt)

    





    if(is.null(encode_track)){ encode_track <- track.gencode()} else{print("Not downloading gencode annotation track")}

    #encode_track$y.field <-
    out_track <- c(encode_track, out_track)

#    out_track@formatting[,1] <- NULL
    out_track$height <- 16
    plot_ranges <- gr.sub2(plot_ranges)
##    print(plot_ranges)

#    print(lapply(out_track, seqlevels))
#    print(out_track)



if(is.null(links)){
    plot.gTrack(c(out_track),windows=plot_ranges,y.heights=y_height)} else{
        values(links)$col ="gray50"
        values(links)$lwd=1
        values(links)$lty=1
        loops <- unique(gr.sub2(links))
                                                  
                                                  
                                                  ##print(loops)
                                                  ##strand(loops) <- 
                                                  ##print(class(loops))
                                                     plot.gTrack(c(out_track),windows=plot_ranges,links=loops,links.feat=values(loops),y.heights=y_height)
                                                 }
   

#    for(i in 1:length(all_gt)){
#        sub <- all_gt[[i]]
                                        #        assign(paste0("gt_",i),sub)


}




collapse_gtrack_list <- function(gtl){
    out_track <- gtl[[1]]
    for(i in 2:length(gtl)){out_track <- c(out_track, gtl[[i]])}
    return(out_track)}
collapse_granges_list <- function(grl){
    out_gr <- grl[[1]]
    out_gr$gr_name <- names(grl)[1]
    for(i in 2:length(grl)){int <- grl[[i]]; int$gr_name <- names(grl)[i]; out_gr <- c(out_gr, int)}
    return(out_gr)}

    
rrbind <- function (..., union = TRUE, as.data.table = FALSE)
    {
        dfs = list(...)
        dfs = dfs[!sapply(dfs, is.null)]
        dfs = dfs[sapply(dfs, ncol) > 0]

        if (any(mix <- sapply(dfs, class) == 'matrix'))
            dfs[mix] = lapply(dfs, as.data.frame)

        names.list = lapply(dfs, names)
        cols = unique(unlist(names.list))
        unshared = lapply(names.list, function(x) setdiff(cols, x))
        ix = which(sapply(dfs, nrow) > 0)
        if (any(sapply(unshared, length) != 0))
            expanded.dts <- lapply(ix, function(x) {
                                       tmp = dfs[[x]]
                                       if (is.data.table(dfs[[x]]))
                                           tmp = as.data.frame(tmp)
                                       tmp[, unshared[[x]]] = NA
                                       return(data.table::as.data.table(as.data.frame(tmp[,
                                                                                          cols])))
                                   })
        else expanded.dts <- lapply(dfs, function(x) as.data.table(as.data.frame(x)[, cols]))

        rout <- tryCatch(rbindlist(expanded.dts), error = function(e) NULL)
        if (is.null(rout))
            rout = data.table::as.data.table(do.call("rbind", lapply(expanded.dts,
                as.data.frame)))
        if (!as.data.table)
            rout = as.data.frame(rout)
        if (!union) {
            shared = setdiff(cols, unique(unlist(unshared)))
            rout = rout[, shared]
        }
        return(rout)
    }


gr.sub2 <- function(gr){
    require(GenomicRanges)
    seqlevels(gr) <- gsub("chr","",seqlevels(gr))
    seqlevels(gr) <- gsub("MT","M",seqlevels(gr))

    return(gr)}

cumulative_link_eval <- function(interesting_links,P_G_link_mat, expression_df,disease_metadata,metadata_col=2,dictionary=NULL) {
    require(reshape2)
    require(dplyr)
    check <- interesting_links %in% rownames(P_G_link_mat)
    if(all(check != TRUE)) { print("Some of the links provided are not in the matrix")} else{
           sub_mat <- P_G_link_mat[interesting_links,]

       }

    melted_mat <- reshape2::melt(as.matrix(sub_mat))
##    print(head(melted_mat))
    melted_mat[,1] <- as.character(melted_mat[,1])
    melted_mat[,2] <- as.character(melted_mat[,2])
    gene <- unique(unlist(lapply(melted_mat[,1], function(x) unlist(strsplit(x,"~"))[2])))
    peaks <- unique(unlist(lapply(melted_mat[,1], function(x) unlist(strsplit(x,"~"))[1])))
    str(gene)
    str(peaks)
    melted_mat$Peaks <- peaks
    melted_mat$Gene <- gene

    if(is.null(dictionary) == FALSE){
        melted_mat$Peaks <- unlist(dictionary[melted_mat$Peaks])} else{print("Doing no translation")}


    melted_mat$Disease <- disease_metadata[as.character(melted_mat[,2]),metadata_col]
    melted_mat$Expression <- unlist(expression_df[gene,as.character(melted_mat$Var2)])
    melted_mat <- group_by(melted_mat,Var2) %>% mutate(count=sum(value))
    melted_mat2 <- melted_mat %>% dplyr::filter(value>0) %>% mutate(Links=paste0(sort(unique(Peaks)), collapse="|")) %>% data.frame
    return(melted_mat2)
}
make_expression_metric_cor_plot <- function(expression_df, metric_df,metric_df_name=NULL ,cluster_df,metric_df2=NULL,metric_df2_name=NULL){
##For making correlation df of super TF outdegree with expression
    require(reshape2)
    if(is.null(metric_df2)){
    common <- intersect(rownames(expression_df), rownames(metric_df))
    sub_exp <- melt(as.matrix(expression_df[common,]))
    sub_met <- melt(as.matrix(metric_df[common,]),id.vars=c("TF","Sample"),variable.name="Metric")

    colnames(sub_exp) <- c("TF","Sample","Expression")
    if(is.null(metric_df_name)){
        colnames(sub_met) <- c("TF","Sample","Metric")} else {
            colnames(sub_met) <- c("TF","Sample",metric_df_name)}

    sub_exp$combo <- paste0(as.character(sub_exp[,1]),"-",as.character(sub_exp[,2]))
    sub_met$combo <- paste0(as.character(sub_met[,1]),"-",as.character(sub_met[,2]))


    out <- merge(sub_exp,sub_met, by="combo")

} else {


        common <- Reduce(intersect,list(rownames(expression_df), rownames(metric_df),rownames(metric_df2)))
        sub_exp <- melt(as.matrix(expression_df[common,]))
        sub_met <- melt(as.matrix(metric_df[common,]),id.vars=c("TF","Sample"),variable.name="Metric")
        sub_met2 <- melt(as.matrix(metric_df2[common,]),id.vars=c("TF","Sample"),variable.name="Metric2")

        colnames(sub_exp) <- c("TF","Sample","Expression")
        if(is.null(metric_df_name)){
        colnames(sub_met) <- c("TF","Sample","Metric")} else {
            colnames(sub_met) <- c("TF","Sample",metric_df_name)}
        if(is.null(metric_df2_name)){
        colnames(sub_met2) <- c("TF","Sample","Metric2")} else {
            colnames(sub_met2) <- c("TF","Sample",metric_df2_name)}

        sub_exp$combo <- paste0(as.character(sub_exp[,1]),"-",as.character(sub_exp[,2]))
        sub_met$combo <- paste0(as.character(sub_met[,1]),"-",as.character(sub_met[,2]))
        sub_met2$combo <- paste0(as.character(sub_met2[,1]),"-",as.character(sub_met2[,2]))

        int <- merge(sub_exp,sub_met, by="combo")
        out <- merge(int,sub_met2,by="combo")



    }


        out <- out[,c(1,2,3,4,7,10)]
    out[,2] <- as.character(out[,2])
    out[,3] <- as.character(out[,3])
    colnames(out)[1:4] <- c("Combo","TF","Sample","Expression")

    out$Super_TF <- NA
    out$TCGA <- NA
    out$Cluster <- NA
    for(i in unique(out$Sample)){
        index <- which(out$Sample == i)

        TF <- out[index,"TF"]
        signif <- unlist(strsplit(cluster_df[i,"Signif_drivers"],"-"))


        bool <- TF %in% signif

        out[index,"Super_TF"] <- bool
        out[index,"TCGA"] <- cluster_df[i,"Subtype"]
        out[index,"Cluster"] <- cluster_df[i,"Cluster"]


        write.table(table(out$Super_TF),"test_cor.txt")
    }
    out$Size <- 0.05+as.numeric(out$Super_TF)
    out$font <- out$Size/4
    return(out)}



plot_expression_metric_cor_heatmap <- function(metric_expression_df,TF,by_cluster=FALSE,rank=TRUE,highlight=FALSE){
    require(dplyr)
    require(reshape2)
    require(ggplot2)

    sub_mat <- metric_expression_df[,2:9]
    sub_mat <- sub_mat[which(sub_mat$TF == TF),]

    zscore_func <- function(vec){
        out <- (vec-mean(vec,na.rm=TRUE))/sd(vec,na.rm=TRUE)
        return(out)}
                                        #    str(sub_mat)
    if(rank == TRUE){
        sub_mat <- sub_mat %>% mutate_if(is.numeric, funs(rank))} else{ sub_mat <- sub_mat %>% mutate_if(is.numeric, funs(zscore_func))}



#    sub_mat$Expression <- rank(sub_mat$Expression)
#    sub_mat$Outdegree <- rank(sub_mat$Outdegree)
##    sub_mat$Outdegree <- rank(sub_mat$Outdegree)

    class_id <- unlist(lapply(1:ncol(sub_mat), function(x) class(sub_mat[,x])))
    index <- which(class_id == "numeric")
    sub_mat$TF_Score <- apply(sub_mat[,index],1, mean)

#    str(sub_mat)

    order_x <- hclust(dist(as.matrix(sub_mat[,c("TF_Score")])))$order


    sub_mat <- sub_mat[order_x,]




    melted_mat <- melt(sub_mat)
#    str(melted_mat)

                                        #    for(i in 1:2){ melted_mat[,i] <- as.character(melted_mat[,i])}
    #melted_mat <- filter(melted_mat,TF==TF)
    melted_mat$Sample <- factor(melted_mat$Sample, levels=melted_mat$Sample[order(melted_mat$value[melted_mat$variable=="TF_Score"])])
    if(highlight==FALSE){
        p <- ggplot(melted_mat, aes(x=Sample, y=variable, fill=value))+geom_raster()+scale_fill_distiller(palette = "Spectral")+labs(title=paste0("Comparison of Metrics for ",TF," across samples"))+theme(axis.text.x=element_text(face='bold',size=4),axis.text.y=element_text(face='bold',size=15),plot.title = element_text(size = rel(2)),strip.text.x = element_text(colour = "black", face = "bold",size=16))} else if (highlight==TRUE){ print("Not working yet")}

    p2 <- p+facet_grid(.~TCGA,scales="free")

    print(p)
    print(p2)
    if(by_cluster==TRUE){

        melted_mat2 <- melted_mat %>% group_by(Cluster) %>% mutate(Cluster2=paste0(unique(TCGA),collapse="\n")) %>% data.frame
        melted_mat2$Cluster <- paste0(melted_mat2$Cluster,"| ",melted_mat2$Cluster2)
        melted_mat2$Sample <- factor(melted_mat2$Sample, levels=melted_mat2$Sample[order(melted_mat2$value[melted_mat2$variable=="TF_Score"])])

        out <- data.frame()

#        for(i in unique(melted_mat2$Cluster)){ sub <- melted_mat2[which(melted_mat2$Cluster == i),]; cast(


                                                                            p3 <- ggplot(melted_mat2, aes(x=Sample, y=variable, fill=value))+geom_raster()+scale_fill_distiller(palette = "Spectral")+labs(title=paste0("Comparison of Metrics for ",TF," across samples"))+theme(axis.text.x=element_text(face='bold',size=4),axis.text.y=element_text(face='bold',size=15),plot.title = element_text(size = rel(2)),strip.text.x = element_text(colour = "black", face = "bold",size=14,hjust=0))+facet_grid(.~Cluster,scales="free"); print(p3)}

#    return(melted_mat)
}


normalize_P_G_mat <- function(Peak_gene_link_mat,peak_height_mat, cutoff=0, quantile=NULL){
    mat <- Peak_gene_link_mat
    peaks <- unlist(lapply(rownames(mat), function(x) unlist(strsplit(x,"~"))[1]))

    if(is.null(quantile) == TRUE){
        link_count <- rowSums(mat)
        peak_count <- rowSums(peak_height_mat[peaks,]>cutoff)} else{





            peak_list <- get_real_peaks(peak_height_mat,quantile)

            true_peak_list <- lapply(peak_list, function(x) intersect(x,peaks))

            peak_count <- as.data.frame(table(unlist(true_peak_list)))

            rownames(peak_count) <- as.character(peak_count[,1])

            peak_count[,1] <- NULL
            name_id <- rownames(peak_count)
            peak_count <- unlist(peak_count)
            names(peak_count) <- name_id


            mat2 <- mat


            for(i in names(true_peak_list)){

                index <- which(peaks %in% true_peak_list[[i]])
                index2 <- which(mat2[,i] == 1)
                final_index <- intersect(index,index2)
                mat2[,i] <- 0
                mat2[final_index,i] <- 1
                link_count <- rowSums(mat2)

            }}


    link_freq <- link_count/peak_count[peaks]

    return(link_freq)}


get_metric_cluster_composition <- function(metric_tsne_df){
    require(reshape2)
    sub <- metric_tsne_df
    test_group <- do.call(rbind,lapply(unique(sub$Cluster), function(x) cbind(paste0("Cluster-",x),as.data.frame(table(sub[which(sub$Cluster == x),"Subtype"]),stringsAsFactors=F)))); colnames(test_group) <- c("Cluster","Disease","Count")
    test_group[,1] <- as.character(test_group[,1])


    out <- data.frame()
    for(i in unique(as.character(sub$Cluster))){
        index <- which(test_group$Cluster == paste0("Cluster-",i))

        int <- test_group[index,]
        supers <- unique(unlist(strsplit(sub[which(sub$Cluster==i),"Cluster_drivers"],"-")))


        int <- do.call(rbind,lapply(supers,function(x) rbind(cbind(int,x), stringsAsFactors=F)))
        int$TF <- as.character(int$x)
        int$x <- NULL
        out <- rbind(out, int)
        out <- as.data.frame(out, stringsAsFactors=F)


    }
        out <- out %>% group_by(TF) %>% mutate( N_Clusters=length(unique(Cluster)))
    out <- out %>% mutate(N_Disease=length(unique(Disease)))
    out <- out %>% mutate(N_TF=sum(Count))  %>% data.frame
    return(out)}

NKX_correction <- function(dir, extension=".rds"){
    require(igraph)
    file_list <- list.files(dir,extension)
    for(i in file_list){
        if(extension ==".rds"){
        infile <- readRDS(paste0(dir,i))
        el <- get.edgelist(infile)
    } else if (extension == ".txt") { infile <- read.table(paste0(dir,i),sep='\t',stringsAsFactors=F)
                                      el <- infile }
        index_mat <- apply(el, 2, function(x) grep("NKX", x))
        for(j in 1:2){
            index <- index_mat[[j]]
            el[index,j] <- gsub("-","",el[index,j])}
        graph_out <- graph_from_edgelist(as.matrix(el[,c(1,2)]))
        saveRDS(graph_out, paste0(dir,gsub(".txt",".rds",i)))}}

get_optimal_dbscan_eps <- function(tsne_df,min_samples=3,opt_count=1,kmin_samples=3){
    require(dbscan)
    require(splines)
    require(inflection)
    require(dplyr)
    sub_df <- tsne_df[,grep("Dim", colnames(tsne_df))]
    knn_df <- sort(kNNdist(sub_df,kmin_samples))
    knn_df <- data.frame(1:length(knn_df),knn_df)

    diff <- derivative(knn_df[,1], knn_df[,2])
   # diff
    #for(i in 2:nrow(knn_df)){ xdiff <- knn_df[i,1]-knn_df[i-1,1]; ydiff <- knn_df[i,2]-knn_df[i-1,2]; diff <- c(diff, (ydiff/xdiff))}
    #diff[which(is.infinite(diff) == TRUE)] <- 0
                                        #diff2 ==diff
#    knn_df[,2] <- round(knn_df[,2],3)
    knn_df$dv_dt <- diff
    #blacklist <- c(1:5, rev(1:nrow(knn_df))[1:3])
    blacklist <- 1:5
    #str(blacklist)
    knn_df <- knn_df[setdiff(1:nrow(knn_df),blacklist),]
    colnames(knn_df) <- c("Sample","Distance","Diff")
    knn_df$Dist_Bisection <- unlist(lapply(1:nrow(knn_df), function(x) dist2d(c(knn_df[x,1],knn_df[x,2]),c(knn_df[1,1],knn_df[1,2]),c(knn_df[nrow(knn_df),1],knn_df[nrow(knn_df),2]))))
    knn_df$Diff_Bisection <- unlist(lapply(1:nrow(knn_df), function(x) dist2d(c(knn_df[x,1],knn_df[x,3]),c(knn_df[1,1],knn_df[1,3]),c(knn_df[nrow(knn_df),1],knn_df[nrow(knn_df),3]))))
    
    knn_df$Elbow_Diff <- bede(knn_df$Sample, knn_df$Diff,0)$iplast
    knn_df$Elbow_bisection <- knn_df[which(knn_df$Dist_Bisection == max(knn_df$Dist_Bisection)),"Sample"]
    knn_df$Elbow_bisection_diff <- knn_df[which(knn_df$Diff_Bisection == max(knn_df$Diff_Bisection)),"Sample"]
    
    
    optimal_dist <- data.frame(stringsAsFactors=F)
    optimal_dist2 <- data.frame(stringsAsFactors=F)
    for(i in seq(0.3,1,0.05)){
        spline_df <- data.frame(smooth.spline(knn_df$Sample, knn_df$Distance)$x,smooth.spline(knn_df$Sample, knn_df$Distance,spar=i)$y,stringsAsFactors=F)
        colnames(spline_df) <- c("Sample","Distance")
        spline_df$span <- paste0("Span:",i)
        spline_df$Elbow <- bede(spline_df$Sample, spline_df$Distance,0)$iplast 

        optimal_dist <- rbind(optimal_dist,spline_df,stringsAsFactors=F)}
#    for(i in seq(0.3,1,0.05)){
#        spline_df2 <- data.frame(smooth.spline(knn_df$Sample, knn_df$Diff)$x,smooth.spline(knn_df$Sample, knn_df$Diff,spar=i)$y,stringsAsFactors=F)
#        colnames(spline_df2) <- c("Sample","Diff")
#        spline_df2$span <- paste0("Span:",i)
#        spline_df2$Diff2 <- derivative(spline_df2$Sample, spline_df2$Diff)
        
#        spline_df2$Elbow <- bede(spline_df2$Sample, spline_df2$Diff,0)$iplast
#        spline_df2$ElbowD2 <- spline_df2[which(spline_df2$Diff2 == max(spline_df2$Diff2))[1],"Sample"]
        #print(unique(spline_df2$Elbow))
        #print(unique(spline_df2$ElbowD2))
#        if(is.na(spline_df2$Elbow==TRUE)

#        optimal_dist2 <- rbind(optimal_dist2,spline_df2,stringsAsFactors=F)}
    
    optimal_dist$Method <- "Direct"
#    optimal_dist2$Method <- "Derivative"
    #print(table(is.na(optimal_dist2$ElbowD2)))
   #print(unique(optimal_dist2[c("Elbow","ElbowD2")]))
    interesting_opts <- unique(optimal_dist[c("Elbow","Method")])
#    interesting_opts <- rbind(interesting_opts,unique(optimal_dist2[c("Elbow","Method")]))
    interesting_opts <- interesting_opts[which(is.na(interesting_opts$Elbow) !=TRUE),]
    interesting_opts$Elbow <- round(interesting_opts$Elbow)
    knn_df$Method <- interesting_opts[1,"Method"]
    
    #print(unique(optimal_dist2[c("Elbow","span")]))
    
    final_opts <- knn_df[which(knn_df$Sample == knn_df$Elbow_bisection),]
    #str(final_opts)
    final_eps <- round(mean(final_opts$Distance[1:opt_count],na.rm=TRUE),digits=3)
    if(is.na(final_eps) == TRUE){ print("Having difficulty finding an optimal root")
                                  saveRDS(knn_df,"knn_df.rds")
                                  saveRDS(optimal_dist,"optimal_dist.rds")
#                                  saveRDS(optimal_dist2,"optimal_dist2.rds")

                                  #str(interesting_opts)
                                  #str(final_opts)
                                  if(is.na(unique(knn_df$Elbow_bisection)) == FALSE) {
                                      print("Using elbow diff")
                                      rounded <- round(unique(knn_df$Elbow_bisection))
                                      str(dplyr::filter(knn_df, Sample == rounded)$Distance)
                                      knn_df$Elbow <- rounded
                                      print(unique(knn_df$Elbow))
                                      knn_df$Method <- "Direct_HG"
                                      final_eps <- round(dplyr::filter(knn_df, Sample == rounded)$Distance,3)} else { print("Using highest differential")
                                                                                                         final_eps <- knn_df[order(-knn_df$Diff_bisection),"Distance"][1]}
                                  print("Best guess is")
                                  print(final_eps)
                              } else{ print(final_eps)}
                                        #str(knn_df[order(-knn_df$Diff),])
    #print(paste0("D2UIK gives ", filter(knn_df, Sample==d2uik(knn_df$Sample, knn_df$Distance))$Distance))
   # knn_df$Elbow_d2uik <- d2uik(knn_df$Sample, knn_df$Distance)
   # knn_df$Distance_d2uik <- filter(knn_df, Sample==d2uik(knn_df$Sample, knn_df$Distance))$Distance
    saveRDS(knn_df,"knn_df.rds")
    saveRDS(optimal_dist,"optimal_dist.rds")
   # saveRDS(optimal_dist2,"optimal_dist2.rds")

    return(final_eps)}


get_optimal_k_clusters <- function(tsne_df, min_k=2,max_k=40, method=c("kmeans","GMM","kmedoids","hclust")){
    require(cluster)
    require(mclust)
    max_k <- min(c(max_k,nrow(tsne_df)-1))
    sub_df <- tsne_df[,grep("Dim", colnames(tsne_df))]

    if(method=="kmeans"){

        dis <- dist(sub_df)

    out_df <- data.frame(stringsAsFactors=F)
    for(i in min_k:max_k){
        kmeans_df <- kmeans(sub_df,i)

        wss <- kmeans_df$tot.withinss
        sil <- silhouette(kmeans_df$cluster,dis)

        out_df <- rbind(out_df,cbind(i,mean(sil[,3]),wss))
    }
#        str(out_df)
    colnames(out_df) <- c("Cluster_Count", "Sil_Size","WithinSS")
    diff <- 0; for(i in 2:nrow(out_df)){ ydiff <- out_df[i,3]-out_df[i-1,3]; xdiff <- out_df[i,1]-out_df[i-1,1]; diff <- c(diff, (ydiff/xdiff))}
        out_df$Elbow_optimization <- diff
#        print(out_df)
        sil_opt <- out_df[order(-out_df$Sil_Size),1][1:5]
#        str(out_df[sil_opt,])
#        str(sil_opt)
        elbow_opt <- out_df[order(-out_df$Elbow_optimization),1][1:5]
#        str(out_df[elbow_opt,])
#        print(out_df[order(-out_df$Elbow_optimization),][1:10,])
#        print(out_df[order(-out_df$Sil_Size),][1:10,])
        optimal_k <- max(intersect(sil_opt,elbow_opt))
#        str(optimal_k)
        if(is.finite(optimal_k) ==FALSE){ optimal_k <- sil_opt[1]
                                        #            print(paste0("No concordance between elbow and silhouette methods"))
                                      }
    }
    else if(method == "GMM"){
        out_df <- data.frame(stringsAsFactors=F)
        for(i in min_k:max_k){
            gmm_df <- Mclust(sub_df,i,verbose=F)
            gmm_k <- log(-1*gmm_df$bic)
            out_df <- rbind(out_df,cbind(i, gmm_k))
        }
        colnames(out_df) <- c("Cluster_Count","BIC")
        gmm_opt <- out_df[order(out_df$BIC),1][1:5]
        optimal_k <- gmm_opt[1]
        sil_opt <- gmm_opt[1]}
    else if(method=="kmedoids"){
                dis <- dist(sub_df)

    out_df <- data.frame(stringsAsFactors=F)
    for(i in min_k:max_k){
        sil <- pam(sub_df,k=i,diss=FALSE)$silinfo$avg.width
        out_df <- rbind(out_df, cbind(i,sil))
    }
                colnames(out_df) <- c("Cluster_Count", "Sil_Size")
                sil_opt <- out_df[order(-out_df$Sil_Size),1][1:5]
                optimal_k <- sil_opt[1]
            }
    else if(method =="hclust"){
        out_df <- data.frame(stringsAsFactors=F)
        for(i in min_k:max_k){
            dis <- dist(sub_df)
            sub_df$Cluster <- cutree(hclust(dist(sub_df)),k=i)
            sil <- mean(silhouette(sub_df$Cluster,dist(sub_df))[,3])
            out_df <- rbind(out_df, cbind(i,sil))}
            colnames(out_df) <- c("Cluster_Count", "Sil_Size")

            sil_opt <- out_df[order(-out_df$Sil_Size),1][1:5]
            optimal_k <- sil_opt[1]
    }



    output <- list(optimal_k, out_df,sil_opt[1])
    names(output) <- c("Optimal_K","Optimization_df","Optimal_K_Silhouette")
    return(output)}

compare_expression_to_random <- function(expression_df,gene_list,gene_list_label="Test_Set",metadata_file=NULL,metadata_column=2){
    require(reshape2)
    check <- setdiff(gene_list,rownames(expression_df))
    if(length(check >=1)){
           print(paste0(paste0(check,collapse=NULL)," missing from expression dataset"))
           sub_exp <- expression_df[setdiff(gene_list,check),]
       } else{
    sub_exp <- expression_df[gene_list,]}
    sub_random <- expression_df[sample(setdiff(rownames(expression_df),gene_list),length(gene_list)*3),]
    melt_exp <- reshape2::melt(as.matrix(sub_exp))
    melt_random <- reshape2::melt(as.matrix(sub_random))
    melt_exp$Status <- gene_list_label
    melt_random$Status <- "Random"

    if(is.null(metadata_file)){
        out_df <- rbind(melt_exp,melt_random)
    } else{
        melt_exp$Metadata <- metadata_file[as.character(melt_exp$Var2),metadata_column]
        melt_random$Metadata <- metadata_file[as.character(melt_random$Var2),metadata_column]
        out_df <- rbind(melt_exp,melt_random)}

    out_df[is.na(out_df)] <- 0
    return(out_df)}

get_cluster_ids_from_metric_list <- function(metric_list,nested=TRUE){
    if(nested == TRUE){ sub_list <- lapply(metric_list, function(x) x[[2]])}
    else if(nested==FALSE){ sub_list <- metric_list}

    common <- Reduce(intersect, lapply(sub_list, function(x) rownames(x)))
    str(common)
    sub_list <- lapply(sub_list, function(x) x[common,])

    cluster_ids <- do.call(cbind,lapply(sub_list, function(x) x["Cluster"]))
    colnames(cluster_ids) <- names(metric_list)
    rownames(cluster_ids) <- common
    return(cluster_ids)}

Calculate_cluster_similarity <- function(cluster_id_df, metadata_df, column=2,method=c("simple","complex","homogen"),both_ways=TRUE){
    df <- cluster_id_df


    comparisons <- combn(colnames(df), 2,simplify=FALSE)
    if(both_ways ==TRUE){
    comparisons <- c(comparisons, lapply(comparisons, rev))} else{ comparisons <- comparisons}


    final_df <- data.frame()

    for(k in comparisons){

        sub_df <- df[k]
        colnames(sub_df) <- c("Metric1","Metric2")
        sub_df$Metric <- paste0(colnames(df[k]),collapse="~")
        sub_df$Agreement <- 0
        Method <- unique(sub_df$Metric)
        print(Method)

        sub_df$Metadata <- metadata_df[rownames(sub_df),column]
        if(method=="simple"){
        for(i in 1:nrow(sub_df)){

            cohort1 <- sub_df[which(sub_df[,1] == sub_df[i,1]),]
            agreement <- length(unique(cohort1[,1]))/length(unique(cohort1[,2]))
            sub_df[i,"Agreement"] <- agreement

        }
        final_df <- rbind(final_df, sub_df, stringsAsFactors=F)
    } else if(method =="complex"){
            for(i in 1:nrow(sub_df)){

                cohort1 <- sub_df[which(sub_df[,1] == sub_df[i,1]),]
                cohort2 <- sub_df[which(sub_df[,2] %in% unique(cohort1[,2])),]

            agreement <- length(cohort1[,1])/length(cohort2[,2])
            sub_df[i,"Agreement"] <- agreement
            }
            final_df <- rbind(final_df, sub_df, stringsAsFactors=F)

        } else if(method=="homogen"){
                sub_df <- df[k]
                str(sub_df)
                str(rownames(sub_df))
                out <- data.frame(stringsAsFactors=F)
                for(i in unique(sub_df[,1])){
                    int1 <- sub_df[which(sub_df[,1] ==i),]
                    secondary <- unique(int1[,2])
                    int2 <- sub_df[which(sub_df[,2] %in% secondary),]
                    agreement <- sum(as.numeric(rownames(int1) %in% rownames(int2)))/nrow(int2)
                    int1$Agreement <- agreement
                    int1$Method <- Method

                    out <- rbind(out, int1)
                   

                }
                final_df <- rbind(final_df,out)
        
            }}

    return(final_df)}


get_cluster_specific_edges <- function(sample_type, metadata_file="/pbtech_mounts/homes024/anf2034/TCGA_Pan_Can_metadata.csv",metadata_column =2, TF_list="/pbtech_mounts/homes024/anf2034/TF_list_Erica_networks.txt",matrix_dir="/pbtech_mounts/homes024/anf2034/ATAC_Seq_Project/Target_Matrices/",num_rows=NULL,min_samples_cutoff=2,limit_to_metadata=TRUE){
    source("/home/anf2034/Andre_F_functions.R")
    if( is.na(file.info(TF_list)$size)==TRUE) { TF_list <- TF_list } else if (is.na(file.info(TF_list)$size)==FALSE){ TF_list <- readLines(TF_list)}
    str(TF_list)
    if(is.character(metadata_file)){
        meta <- read.table(metadata_file,header=T,sep=',',stringsAsFactors=F)} else if (is.data.frame(metadata_file)){ meta <- metadata_file}
#    str(metadata_file)
   print(sample_type)
    in_samples <- meta[which(meta[metadata_column] == sample_type),"Sample"]
#    str(in_samples)
    if(length(in_samples) <=1){ print("This cohort is too small for reliable analysis")} else{

                                        #    subtype_df <- data.frame(stringsAsFactors=F)
        ratio_list <- list()
    for(i in TF_list){
        in_mat <- readRDS(paste0(matrix_dir,i,"_target_matrix.rds"))
        if(limit_to_metadata==TRUE){
            in_mat <- in_mat[intersect(rownames(in_mat), meta$Sample),]} else {in_mat <- in_mat}
        in_samples_real <- intersect(in_samples,rownames(in_mat))
#        str(in_samples_real)
        sub_mat <- in_mat[in_samples_real,]
        in_mat <- in_mat[setdiff(rownames(in_mat),in_samples_real),]
        inter_summary <- round(colSums(in_mat)/nrow(in_mat),digits=3)
        intra_summary <- round(colSums(sub_mat)/nrow(sub_mat),digits=3)
        ratio <- intra_summary/inter_summary
        ratio <- data.frame(inter_summary,intra_summary,ratio,stringsAsFactors=F)

        ratio$Target <- rownames(ratio)
        ratio$TF <- i
        ratio$Disease <- sample_type
        colnames(ratio)[3] <- "Ratio"
        ratio$Ratio <- as.numeric(gsub(Inf,10000,ratio$Ratio))
#        str(ratio)
        ratio <- ratio[which(ratio$intra_summary >= (min_samples_cutoff/nrow(sub_mat))),]
        if(is.null(num_rows)) { ratio <- ratio} else{ ratio <- ratio[order(-ratio$Ratio),]
                                                      ratio <- ratio[1:num_rows,]}
        ratio_list[[i]] <- ratio }

        return(ratio_list)}




#        inter_summary <- 100*table(inter_summary>cutoff)/sum(table(inter_summary>cutoff))
#        print(inter_summary)
#        intra_summary <- 100*table(intra_summary>cutoff)/sum(table(intra_summary>cutoff))
#        print(intra_summary)

#        subtype_df <- rbind(subtype_df,cbind(inter_summary["TRUE"],intra_summary["TRUE"],sample_type,i),stringsAsFactors=F)
#        print(paste0("Done ",i))
}

#    colnames(subtype_df) <- c("Inter_type","Intra_type","Disease_Type","TF")
#    saveRDS(subtype_df, paste0(outdir,sample_type,"_target_variability.rds"))
#    write.table(subtype_df, paste0(outdir,sample_type,"_target_variability.txt"),quote=F, sep='\t')
#    if(return_df == TRUE){return(subtype_df)} else{ print(paste0("Finished ",sample_type))}}






rewiring_auc_clusters <- function(dir,tsne_df,tsne_column="Disease"){
    file_list <- list.files(dir,"_matrix.rds")
#    file_list <- file_list[15:18]


    out_list <- list()
    for(i in file_list){
        print(i)
        tsne_df$Check <- NA
        tsne_df$Check <- unlist(lapply(tsne_df$Cluster_drivers, function(x) i %in% unlist(strsplit(x,"-"))))
#        str(tsne_df)
        in_mat <- readRDS(paste0(dir,i))
#        str(in_mat)
        index_list <- lapply(unique(tsne_df[,tsne_column]), function(x) tsne_df[which(tsne_df[,tsne_column] == x),c("Sample","Check")])
#        str(index_list)
        names(index_list) <- unique(tsne_df[,tsne_column])
        cluster_mat <- data.frame(0)
        for(j in 1:length(index_list)){
            index <- index_list[[j]][,1]
            sub_mat <- as.data.frame(in_mat[index,])
#            str(sub_mat)
            if(ncol(sub_mat) >1){

                targets <- colnames(sub_mat)[which(colSums(sub_mat) >1)]
                if (length(targets) >1){
                    target_df <- data.frame(1,targets,stringsAsFactors=F)
                    rownames(target_df) <- targets
                    cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row")
                } else{ targets <- rownames(cluster_mat)
                        target_df <- data.frame(0,targets,stringsAsFactors=F)
                        rownames(target_df) <- targets

                        cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row")
                    }

            target_df$targets <- NULL

#                str(rownames(cluster_mat))


        } else{

            targets <- rownames(sub_mat)[which(sub_mat >=1)]
            if(length(targets) >1){
                target_df <- data.frame(1,targets,stringsAsFactors=F)} else{
                    targets <- rownames(cluster_mat)
                    target_df <- data.frame(0,targets,stringsAsFactors=F)}
            rownames(target_df) <- targets
            #str(target_df)
            target_df$targets <- NULL
            cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row") }}
        cluster_mat <- cluster_mat[-1,-1]

        colnames(cluster_mat) <- names(index_list)
        saveRDS(cluster_mat,paste0(dir,gsub(".rds","_clusters.rds",i)))
        cluster_mat <- t(cluster_mat)
        cluster_mat[is.na(cluster_mat)] <- 0
#        print(table(colSums(cluster_mat)<10))



        initial <- names(sort(rowSums(cluster_mat)))[1]
        blacklist <- c(initial)
        target_len <- names(which(cluster_mat[initial,] == 1))
        target_all_len <- length(target_len)

        for(j in 2:nrow(cluster_mat)){
            sub_mat <- cluster_mat[setdiff(rownames(cluster_mat),blacklist),]
            if(is.null(nrow(sub_mat))){target_index <- setdiff(rownames(cluster_mat),blacklist) } else{ sub_mat <- sub_mat
            target_list <- apply(sub_mat, 1, function(x) colnames(sub_mat)[which(x == 1)])
            target_list <- sort(unlist(lapply(target_list, function(x) length(setdiff(x, target_len)))))



            target_index <- names(target_list)[1]}

            blacklist <- c(blacklist,target_index)
            target_vec <- names(which(cluster_mat[target_index,] == 1))


            target_len <- union(target_len,target_vec)
            target_all_len <- c(target_all_len,length(target_len))
#            str(blacklist)



        }
        target_all_len <- as.data.frame(target_all_len,stringsAsFactors=F)
        colnames(target_all_len) <- "Union_Targets"

        target_all_len$Num_Samples <- 1:dim(target_all_len)[1]

        out_list <- c(out_list,list(target_all_len))

 #       print(paste0("Finished ",gsub("_target_matrix.rds","",i)))





    }
    names(out_list) <- gsub("_target_matrix.rds","",file_list)
    return(out_list)
    }

rewiring_auc_clusters_specific <- function(dir,tsne_df, extended=FALSE){
    file_list <- list.files(dir,"_matrix.rds")
#    file_list <- file_list[25:30]


    out_list <- list()
    for(i in file_list){
        TF <- unlist(strsplit(i,"_"))[1]
        print(paste0("Working ",TF))
        tsne_df$Check <- NA
        if(extended == FALSE){
            tsne_df$Check <- unlist(lapply(tsne_df$Cluster_drivers, function(x) TF %in% unlist(strsplit(x,"-"))))} else{
                tsne_df$Check <- unlist(lapply(tsne_df$Signif_drivers, function(x) TF %in% unlist(strsplit(x,"-"))))}
#        print(table(tsne_df$Check))

        in_mat <- readRDS(paste0(dir,i))
#        print(dim(in_mat))
        if(any(tsne_df$Check == TRUE)){
            target_index <- tsne_df[which(tsne_df$Check == TRUE),"Sample"]
 #           str(target_index)
            sub_mat <- in_mat[target_index,]
            in_mat_final <- in_mat[,which(colSums(sub_mat) >=(0.5*nrow(sub_mat)))]
            if(ncol(in_mat)<2){ in_mat <- in_mat
                            } else{ in_mat <- in_mat_final}}
         else{ in_mat <- in_mat}
#        print(dim(in_mat))
        index_list <- lapply(unique(tsne_df$Disease), function(x) tsne_df[which(tsne_df$Disease == x),c("Sample","Check")])
        names(index_list) <- unique(tsne_df$Disease)
        cluster_mat <- data.frame(0)
        for(j in 1:length(index_list)){
            index <- index_list[[j]][,1]
            sub_mat <- as.data.frame(in_mat[index,])
#            str(sub_mat)
            if(ncol(sub_mat) >1){

                targets <- colnames(sub_mat)[which(colSums(sub_mat) >1)]
                if (length(targets) >1){
                    target_df <- data.frame(1,targets,stringsAsFactors=F)
                    rownames(target_df) <- targets
                    cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row")
                } else{ targets <- rownames(cluster_mat)
                        target_df <- data.frame(0,targets,stringsAsFactors=F)
                        rownames(target_df) <- targets

                        cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row")
                    }

            target_df$targets <- NULL

#                str(rownames(cluster_mat))


        } else{

            targets <- rownames(sub_mat)[which(sub_mat >=1)]
            if(length(targets) >1){
                target_df <- data.frame(1,targets,stringsAsFactors=F)} else{
                    targets <- rownames(cluster_mat)
                    target_df <- data.frame(0,targets,stringsAsFactors=F)}
            rownames(target_df) <- targets
            #str(target_df)
            target_df$targets <- NULL
            cluster_mat <- vector_smartmerge(cluster_mat,target_df,"row") }}
        cluster_mat <- cluster_mat[-1,-1]

        colnames(cluster_mat) <- names(index_list)
        saveRDS(cluster_mat,paste0(dir,gsub(".rds","_clusters.rds",i)))
        cluster_mat <- t(cluster_mat)
#        print(table(colSums(cluster_mat)<10))



        initial <- names(sort(rowSums(cluster_mat)))[1]
        blacklist <- c(initial)
        target_len <- names(which(cluster_mat[initial,] == 1))
        target_all_len <- length(target_len)

        for(j in 2:nrow(cluster_mat)){
            sub_mat <- cluster_mat[setdiff(rownames(cluster_mat),blacklist),]
            if(is.null(nrow(sub_mat))){target_index <- setdiff(rownames(cluster_mat),blacklist) } else{ sub_mat <- sub_mat
            target_list <- apply(sub_mat, 1, function(x) colnames(sub_mat)[which(x == 1)])
            target_list <- sort(unlist(lapply(target_list, function(x) length(setdiff(x, target_len)))))



            target_index <- names(target_list)[1]}

            blacklist <- c(blacklist,target_index)
            target_vec <- names(which(cluster_mat[target_index,] == 1))


            target_len <- union(target_len,target_vec)
            target_all_len <- c(target_all_len,length(target_len))
#            str(blacklist)



        }
        target_all_len <- as.data.frame(target_all_len,stringsAsFactors=F)
        colnames(target_all_len) <- "Union_Targets"

        target_all_len$Num_Samples <- 1:dim(target_all_len)[1]

        out_list <- c(out_list,list(target_all_len))

 #       print(paste0("Finished ",gsub("_target_matrix.rds","",i)))





    }
    names(out_list) <- gsub("_target_matrix.rds","",file_list)
    return(out_list)
    }
make_rewiring_df <- function(rewiring_output_list,Super_TF_list){

    for(i in 1:length(rewiring_output_list)){
        sub <- rewiring_output_list[[i]]
        sub$TF <- names(rewiring_output_list)[i]; rewiring_output_list[[i]] <- sub}
    final_output_list <- do.call("rbind", rewiring_output_list)
    final_output_list <- final_output_list %>% group_by(TF) %>% mutate(Norm_Targets=Union_Targets/max(Union_Targets)) %>% data.frame
    final_output_list <- final_output_list %>% group_by(TF) %>% mutate(Norm_Samples=Num_Samples/max(Num_Samples)) %>% data.frame
    final_output_list <- final_output_list %>% group_by(TF) %>% mutate(Area=trapz(Norm_Samples,Norm_Targets)) %>% data.frame
    final_output_list <- final_output_list %>% group_by(TF) %>% mutate(Max_Union_Targets=max(Union_Targets)) %>% data.frame
    final_output_list$Super_TF <- final_output_list$TF %in% Super_TF_list
    final_output_list$Super_TF <- gsub(TRUE,"Yes",final_output_list$Super_TF)
    final_output_list$Super_TF <- gsub(FALSE,"No",final_output_list$Super_TF);final_output_list$Method <- "Cluster_Targets"
    final_output_list$Log_Targets <- log(final_output_list$Max_Union_Targets)
    return(final_output_list)}


calc_metric_stability <- function(indir, outfile,fraction,metric=c("outdegree","indegree","betweenness"),verbose=FALSE){
    require(igraph)

    file_list <- list.files(indir,"_graph.rds")

    metric_df <- data.frame(0)
    rownames(metric_df) <- "Empty"
    str(file_list)

    print(paste0("Calculating ",metric))

    for(i in file_list){
#        str(i)
        graph <- readRDS(paste0(indir,i))
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
    colnames(metric_df) <- gsub("_Final|_PIQ|_graph.rds","",colnames(metric_df))
    metric_df <- metric_df[-1,-1]

saveRDS(metric_df,outfile)
}


run_umap <- function(df,sample_margin="column",config=umap.defaults,seed=20,n_neighbors=3,min_dist=0.1,spread=1){
    require(umap)
    if(sample_margin == "column"){ df <- as.data.frame(t(df),stringsAsFactors=F)
                                   rownames(df) <- gsub("\\.","-",rownames(df))
                               }
    else if (sample_margin == "row") { df <- df}
    config$random_state=seed
    config$n_neighbors=n_neighbors
    config$min_dist <- min_dist
    config$spread <- spread

    out <- umap(df,config)$layout
    out <- as.data.frame(out,stringsAsFactors=F)
    out$Sample <- rownames(df)
    rownames(out) <- out[,3]
    colnames(out) <- c("Dim1","Dim2","Sample")
    return(out)
}

make_clusters_from_metric_list_umap <- function(metric_list,metadata_file,method="dbscan"){
    metric_clusters <- list()
    for(i in 1:length(metric_list)){
        print(paste0("Working ",names(metric_list)[i]))
        df <- metric_list[[i]]

        tsne_out <- run_umap(df)
        print("UMAP done")

        tsne_out <- get_clusters_tsne(tsne_out,method)

        tsne_out <- get_metadata_tsne(tsne_out, metadata_file)
        str(rownames(tsne_out))

        #print(head(tsne_out))
        tsne_out <- get_drivers_from_clusters(df,tsne_out)
        #str(tsne_out)
        metric_clusters <- c(metric_clusters, list(tsne_out)); print(paste0("Finished metric ",i))}
    names(metric_clusters) <- names(metric_list)
    return(metric_clusters)
}

make_clusters_from_metric_list_umap_mc <- function(metric_list,metadata_file,method="dbscan",n_cores=4){
    require(foreach)
    require(doMC)
    registerDoMC(cores=n_cores)

    metric_clusters <- list()
    foreach(i=1:length(metric_list),.combine="c",.multicombine=TRUE) %dopar% {
        print(paste0("Working ",names(metric_list)[i]))
        df <- metric_list[[i]]

        tsne_out <- run_umap(df)
        print("UMAP done")

        tsne_out <- get_clusters_tsne(tsne_out,method)

        tsne_out <- get_metadata_tsne(tsne_out, metadata_file)
        str(rownames(tsne_out))

        #print(head(tsne_out))
        tsne_out <- get_drivers_from_clusters(df,tsne_out)
        #str(tsne_out)
        metric_clusters <- c(metric_clusters, list(tsne_out)); print(paste0("Finished metric ",i))
        return(metric_clusters)
    }
#    names(metric_clusters) <- names(metric_list)
}

math_ops_on_matrix_list <- function(matrix_list, operation){
    check <- length(unique(lapply(matrix_list, dim)))
    if(check >1){ print("These matrices are not conformable")} else{
        check_row <- length(unique(lapply(matrix_list,rownames)))
        check_col <- length(unique(lapply(matrix_list,colnames)))
        if(any(c(check_row,check_col)) >1){ print("Checking Matrices.....")
                                                  print("Matrices are not ordered identically.... Fixing")

                                        } else{ print("Checking Matrices.....")
                                               print("Done!")}

    nrows <- nrow(matrix_list[[1]])
    ncols <- ncol(matrix_list[[1]])

    row_names <- rownames(matrix_list[[1]])
    col_names <- colnames(matrix_list[[1]])


    int_mat <- do.call("cbind", lapply(matrix_list, unlist))
    int_mat_op <- apply(int_mat,1,operation)

    out_mat <- matrix(int_mat_op,nrows,ncols,dimnames=list(row_names,col_names))}
    return(out_mat)}

fishers_exact_vec <- function(vec,universe,label=NULL,alternative_opt=c("two.sided","greater","less"),log_OR=TRUE){
    vec_categories <- unique(vec)
    vec_len <- length(vec)
    int_list <- lapply(vec_categories,function(x) length(which(vec ==x)))
    names(int_list) <- vec_categories
    row1 <- lapply(int_list, function(x) cbind(x,vec_len-x))

    universe_categories <- unique(universe)
    universe_len <- length(universe)
    universe_list <- lapply(universe_categories,function(x) length(which(universe ==x)))
    names(universe_list) <- universe_categories
    row2 <- lapply(universe_list, function(x) cbind(x,universe_len-x))

    common <- intersect(names(row1), names(row2))
    row2 <- lapply(common, function(x) row2[[x]]-row1[[x]])
    names(row2) <- common
    out <- lapply(common, function(x) rbind(row1[[x]],row2[[x]]))
    names(out) <- common

    print(out)
    pvals <- unlist(lapply(out, function(x) fisher.test(x,alternative=alternative_opt)$p.val))
    OR <- unlist(lapply(out, function(x) fisher.test(x,alternative=alternative_opt)$estimate))
    OR_interval <- do.call("rbind",lapply(out, function(x) cbind(fisher.test(x,alternative=alternative_opt)$conf.int[1],fisher.test(x,alternative=alternative_opt)$conf.int[2])))
    ##str(OR_interval)
   

    final <- as.data.frame(cbind(common,pvals,OR,OR_interval,label))
  
##    str(final)
    colnames(final) <- c("Category","P_value","Odds_Ratio","OR_Lower","OR_Upper","Label")
    if(log_OR==TRUE){
        final[2:5] <- apply(final[2:5],2,function(x) as.numeric(x))
        final[3:5] <- apply(final[3:5],2,function(x) log(x))
    final$Is_Log_Odds =TRUE} else{
        final[2:5] <- apply(final[2:5],2, as.numeric)
        final$Is_Log_Odds=FALSE
    }
    return(final)
}

project_onto_tsne <- function(tsne_out, df, margin=c("row","column"), index_df,rank=TRUE,shape=NULL, nudge=-0.5){
    require(ggplot2)
    require(ggrepel)
#    str(tsne_out)

    if(margin == "row"){

        tsne_out[index_df] <- unlist(df[index_df,tsne_out$Sample])} else if (margin == "column"){ tsne_out[index_df] <- unlist(df[tsne_out$Sample,index_df])}

    tsne_out[index_df] <- round(tsne_out[index_df],2)

#    str(tsne_out)
                                        #    print(index_df)
    if(is.null(shape)== TRUE){

        p <- ggplot(tsne_out, aes(Dim1,Dim2))+theme_bw()+geom_point(aes(fill=!!sym(index_df)),shape=21,color="gray50",size=6,show.legend=T)+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))
            ## geom_text_repel(fontface="bold",nudge_y=nudge,size=3.5, aes(color=get(index_df)))
    } else{
        p <- ggplot(tsne_out, aes(Dim1,Dim2))+theme_bw()+geom_point(aes(fill=!!sym(index_df),shape=!!sym(shape)),color="gray50",size=6,show.legend=T)+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+scale_shape_manual(values=c(21:25))
        ## +geom_text_repel(fontface="bold",nudge_y=nudge,size=3.5, aes(color=get(index_df)))
        }

    if(rank==TRUE){ p <- p+scale_fill_gradient(low="red",high="green")} else if (rank==FALSE){
        p <- p+scale_fill_gradient(low="green",high="red")}

    return(p)}


fanying_heatmaps <- function(metric_list,metric_names,cluster,cluster_column,tsne_df,tsne_TFs_column){
    require(ComplexHeatmap)
    require(dplyr)
    require(reshape2)
    require(circlize)


    common_TFs <- Reduce("intersect",lapply(metric_list, rownames))
    common_samples <- Reduce("intersect",lapply(metric_list, colnames))

    metric_list <- lapply(metric_list, function(x) rank_mat(x[common_TFs,common_samples]))
    TF_score <- Reduce("+", metric_list, accumulate=FALSE)/length(metric_list)


    tsne_df <- dplyr::filter(tsne_df, get(cluster_column)==cluster)
    samples <- intersect(tsne_df$Sample,colnames(TF_score))

    TFs <- unlist(strsplit(unique(tsne_df[,tsne_TFs_column]),"-"))
#    str(TFs)
#    str(samples)
    metric_list <- lapply(metric_list, function(x) x[TFs,samples])
    TF_score <- Reduce("+", metric_list, accumulate=FALSE)/length(metric_list)
#    str(metric_list)


#    str(TF_score)


    final_list <- lapply(metric_list, function(x) as.matrix(apply(x,1,mean)))
    TF_score_final <- as.matrix(apply(TF_score,1,median))

#    str(final_list)
#    str(TF_score_final)
    names(final_list) <- metric_names

    col_pal_list <- list(colorRamp2(c(1,340,680),colors=c("blue","white","red"),transparency=0.7),colorRamp2(c(1,340,680),colors=c("green","white","orange"),transparency=0.7),colorRamp2(c(1,340,680),colors=c("yellow","white","purple"),transparency=0.7))
    names(col_pal_list) <- metric_names

    ht_list <- NULL
    for(i in names(final_list)){ ht_list <- ht_list + Heatmap(final_list[[i]],cluster_rows=FALSE,cluster_columns=FALSE,heatmap_legend_param=list(color_bar="continuous"),name=i,col=col_pal_list[[i]],show_row_names=FALSE)}
-    ha <- rowAnnotation(barplot=anno_barplot(TF_score_final,which="row"))
    ht_list <- ht_list+ha
    draw(ht_list)




}

fanying_heatmaps_40 <- function(metric_list,metric_names,cluster,cluster_column,tsne_df,tsne_TFs_column,relative=FALSE,TF_count=10){
    require(ComplexHeatmap)
    require(dplyr)
    require(reshape2)
    require(circlize)


    common_TFs <- Reduce("intersect",lapply(metric_list, rownames))
    common_samples <- Reduce("intersect",lapply(metric_list, colnames))

    metric_list <- lapply(metric_list, function(x) rank_mat(x[common_TFs,common_samples]))
    TF_score <- Reduce("+", metric_list, accumulate=FALSE)/length(metric_list)
    str(TF_score)


    sub_tsne_df <- dplyr::filter(tsne_df, get(cluster_column)==cluster)
    samples <- intersect(sub_tsne_df$Sample,colnames(TF_score))
    if(length(samples) <=1){ print("Not enough samples for this analysis")
                                        #ht_list <- Heatmap(1:10)
                             print(paste0("Working: ",cluster))
                             out_samples <- setdiff(colnames(TF_score), samples)


                             TFs <- unlist(strsplit(unique(sub_tsne_df[,tsne_TFs_column]),"-"))
                             if(length(TFs) < TF_count){ TF_count=length(TFs)} else{TF_count=TF_count}
                             TFs <- TFs[1:TF_count]
                             print(TFs)

                             sub_metric_list <- lapply(metric_list, function(x) x[TFs,samples])
                             out_metric_list <- lapply(metric_list, function(x) x[TFs,out_samples])

                             TF_score <- Reduce("+", sub_metric_list, accumulate=FALSE)/length(sub_metric_list)
                             out_TF_score <- Reduce("+", out_metric_list, accumulate=FALSE)/length(out_metric_list)

                             print("break")

                             if(relative==FALSE){
                                 plot_breaks <- c(1,340,680)
                                 final_list <- sub_metric_list
                                        #str(TF_score)
                                 TF_score_final <- TF_score
                                 ha_bar <- columnAnnotation(TF_Score=anno_barplot(TF_score_final,label=TFs,axis_param = list(gp=gpar(fontsize=20,font=2,direction="reverse"))),height=unit(2.75,"cm"),show_legend=TRUE)
                                        #      print(any(is.na(TF_score_final)))
                                 print(TF_score_final[1:5,])
                                 col_order <- rev(order(TF_score_final))
                                 col_pal_list <- list(colorRamp2(c(1,340,680),colors=c("blue","white","red"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("green","white","orange"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("yellow","white","purple"),transparency=0.3))

                             } else if(relative==TRUE){
                                 plot_breaks <- c(1,340,680)

                                 final_list <- sub_metric_list
                                        #str(TF_score)
                                        #final_out_list <- lapply(out_metric_list, function(x) as.matrix(apply(x,1,mean)))
                                        #final_list <- lapply(1:length(final_list), function(x) final_list[[x]]/final_out_list[[x]])
                                        #TF_score_final <- as.matrix(apply(TF_score,1,mean))
                                 TF_score_final <- cbind(TF_score,as.matrix(apply(out_TF_score,1,mean)))
                                        #print(TF_score_final2[1:5,])
                                 str(TF_score_final)
                                 TF_score_final_int <- t(apply(TF_score_final,1,function(x) x/x[1]))
                                 TF_score_final_int[,1] <- TF_score_final_int[,1]-TF_score_final_int[,2]
                                 TF_score_final_int <- cbind(TF_score_final_int,TF_score_final[,1])

                                 TF_score_final_int[,1] <- TF_score_final_int[,1]*TF_score_final_int[,3]
                                 TF_score_final_int[,2] <- TF_score_final_int[,2]*TF_score_final_int[,3]

                                 TF_score_final <- TF_score_final_int[,c(2,1)]
                                        #print(TF_score_final[1:5,1:2])
                                 final <- apply(TF_score_final,1, sum)
                                 col_order <- rev(order(final))

                                        #TF_score_final <- as.matrix(apply(TF_score,1,mean))
                                 ha_bar <- columnAnnotation(TF_Score=anno_barplot(TF_score_final,gp=gpar(fill=2:3),label=TFs,axis_param = list(gp=gpar(fontsize=20,font=2))),height=unit(2.75,"cm"),show_legend=TRUE)
                                 col_pal_list <- list(colorRamp2(c(1,340,680),colors=c("blue","white","red"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("green","white","orange"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("yellow","white","purple"),transparency=0.3))
                             }






                                        #  str(TF_score_final)
                                        #  print(apply(TF_score_final,2,summary))

                                        #    str(final_list)
                                        #    str(TF_score_final)
                             final_list <- lapply(final_list, function(x) t(x))
                                        #  str(final_list)
                             names(final_list) <- metric_names
                                        #str(TF_score_final)

                                        #  print(TF_score_final[col_order])




        
                             names(col_pal_list) <- metric_names

                             ht_list <- NULL
                             for(i in names(final_list)){ ht_list <- ht_list %v% Heatmap(final_list[[i]],height=unit(3.5,"cm"),row_title = gsub("Ranked_"," ",i),row_title_gp = gpar(font = 2,fontsize=20),row_title_rot=0,rect_gp = gpar(col = "black", lwd = 1.5),column_order= col_order,heatmap_legend_param=list(color_bar="continuous",grid_width=unit(0.75,"cm"),labels_gp = gpar(font = 2),legend_height = unit(2.5, "cm"),at=c(plot_breaks[1],last(plot_breaks)),labels=c("low","high")),name=gsub("_|TFs"," ",i),col=col_pal_list[[i]])}
                                        #ha <- columnAnnotation(TF_Score=anno_barplot(TF_score_final2,gp=gpar(fill=2:3),labels=TFs,height=unit(2.5,"cm")))

                             ha_text <- columnAnnotation(TFs=anno_text(TFs,rot = 30, gp = gpar(fontsize = 18,fontface="bold")))
                             ht_list <- ht_list %v% ha_bar %v% ha_text
                             heatmap_column_names_gp <- TFs
                             ht_list@column_title <- cluster
                             return(ht_list)


                         } else{
                             print(paste0("Working: ",cluster))
                             out_samples <- setdiff(colnames(TF_score), samples)


                             TFs <- unlist(strsplit(unique(sub_tsne_df[,tsne_TFs_column]),"-"))
                             if(length(TFs) < TF_count){ TF_count=length(TFs)} else{TF_count=TF_count}
                             TFs <- TFs[1:TF_count]


                             sub_metric_list <- lapply(metric_list, function(x) x[TFs,samples])
                             out_metric_list <- lapply(metric_list, function(x) x[TFs,out_samples])

                             TF_score <- Reduce("+", sub_metric_list, accumulate=FALSE)/length(sub_metric_list)
                             out_TF_score <- Reduce("+", out_metric_list, accumulate=FALSE)/length(out_metric_list)

                             if(relative==FALSE){
                                 plot_breaks <- c(1,340,680)
                                 final_list <- lapply(sub_metric_list, function(x) as.matrix(apply(x,1,mean)))
                                 TF_score_final <- as.matrix(apply(TF_score,1,mean))
                                 ha_bar <- columnAnnotation(TF_Score=anno_barplot(TF_score_final,label=TFs,axis_param = list(gp=gpar(fontsize=12,font=2,direction="reverse"))),height=unit(2.5,"cm"),show_legend=TRUE)
                                        #      print(any(is.na(TF_score_final)))
                                 print(TF_score_final[1:5,])
                                 col_order <- rev(order(TF_score_final))
                                 col_pal_list <- list(colorRamp2(c(1,340,680),colors=c("blue","white","red"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("green","white","orange"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("yellow","white","purple"),transparency=0.3))

                             } else if(relative==TRUE){
                                 plot_breaks <- c(1,340,680)

                                 final_list <- lapply(sub_metric_list, function(x) as.matrix(apply(x,1,mean)))
                                        #final_out_list <- lapply(out_metric_list, function(x) as.matrix(apply(x,1,mean)))
                                        #final_list <- lapply(1:length(final_list), function(x) final_list[[x]]/final_out_list[[x]])
                                        #TF_score_final <- as.matrix(apply(TF_score,1,mean))
                                 TF_score_final <- cbind(as.matrix(apply(TF_score,1,mean)),as.matrix(apply(out_TF_score,1,mean)))
                                        #print(TF_score_final2[1:5,])
                                 str(TF_score_final)
                                 TF_score_final_int <- t(apply(TF_score_final,1,function(x) x/x[1]))
                                 TF_score_final_int[,1] <- TF_score_final_int[,1]-TF_score_final_int[,2]
                                 TF_score_final_int <- cbind(TF_score_final_int,TF_score_final[,1])

                                 TF_score_final_int[,1] <- TF_score_final_int[,1]*TF_score_final_int[,3]
                                 TF_score_final_int[,2] <- TF_score_final_int[,2]*TF_score_final_int[,3]

                                 TF_score_final <- TF_score_final_int[,c(2,1)]
                                        #print(TF_score_final[1:5,1:2])
                                 final <- apply(TF_score_final,1, sum)
                                 col_order <- rev(order(final))

                                        #TF_score_final <- as.matrix(apply(TF_score,1,mean))
                                 ha_bar <- columnAnnotation(TF_Score=anno_barplot(TF_score_final,gp=gpar(fill=2:3),label=TFs,axis_param = list(gp=gpar(fontsize=12,font=2))),height=unit(2.5,"cm"),show_legend=TRUE)
                                 col_pal_list <- list(colorRamp2(c(1,340,680),colors=c("blue","white","red"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("green","white","orange"),transparency=0.3),colorRamp2(c(1,340,680),colors=c("yellow","white","purple"),transparency=0.3))
                             }






                                        #  str(TF_score_final)
                                        #  print(apply(TF_score_final,2,summary))

                                        #    str(final_list)
                                        #    str(TF_score_final)
                             final_list <- lapply(final_list, function(x) t(x))
                                        #  str(final_list)
                             names(final_list) <- metric_names
                                        #str(TF_score_final)

                                        #  print(TF_score_final[col_order])




                             names(col_pal_list) <- metric_names

                             ht_list <- NULL
                             for(i in names(final_list)){ ht_list <- ht_list %v% Heatmap(final_list[[i]],height=unit(3.5,"cm"),row_title = gsub("Ranked_"," ",i),row_title_gp = gpar(font = 2,fontsize=20),row_title_rot=0,rect_gp = gpar(col = "black", lwd = 1.5),column_order= col_order,heatmap_legend_param=list(color_bar="continuous",grid_width=unit(0.75,"cm"),labels_gp = gpar(font = 2),legend_height = unit(2.5, "cm"),at=c(plot_breaks[1],last(plot_breaks)),labels=c("low","high")),name=gsub("_|TFs"," ",i),col=col_pal_list[[i]])}
                                        #ha <- columnAnnotation(TF_Score=anno_barplot(TF_score_final2,gp=gpar(fill=2:3),labels=TFs,height=unit(2.5,"cm")))

                             ha_text <- columnAnnotation(TFs=anno_text(TFs,rot = 30, gp = gpar(fontsize = 18,fontface="bold")))
                             ht_list <- ht_list %v% ha_bar %v% ha_text
                             heatmap_column_names_gp <- TFs
                             ht_list@column_title <- cluster
                             return(ht_list)




                         }
}



get_coessential_genes <- function(Achilles_df, metadata_df=NULL,subset=NULL, genes_of_interest=rownames(Achilles_df),precomputed_mat=NULL,corr_cutoff=0.95,cor_method=c("pearson","spearman")){
    require(Hmisc)
    require(dplyr)
    require(reshape2)
    print(dim(Achilles_df))
   
    if(is.null(metadata_df) || is.null(subset)){
        print("Working with all data")
        Achilles_df <- Achilles_df} else {
            Achilles_df <- Achilles_df[,intersect(colnames(Achilles_df),metadata_df[,1])]
            metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(Achilles_df)),]
            index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
            print(index)

            if(length(index) <5 & length(index)>0){ print(paste0("Not enough observations for reliable analysis, proceeding with subset + ",5-length(index)," random samples"))
                                  index <- c(index, sample(setdiff(1:nrow(metadata_df), index),5-length(index)))

                                    } else if (length(index) == 0){ print("No cell lines match your search criteria, working with all data")
                                                                    #Achilles_df <- Achilles_df
                                                                } else{
            Achilles_df <- Achilles_df[metadata_df[index,1]]}}
   
    if(is.null(precomputed_mat)){
        print("Calculating correlation matrix")
        cormat <- rcorr(t(Achilles_df),type=cor_method)} else { print("Working with precomputed matrix")
                                                cormat <- precomputed_mat}
    print("Finished Correlation matrix")

    gc()
    cor_vals <- reshape2::melt(cormat[[1]])
    cor_vals[,1:2] <- apply(cor_vals[,1:2],2,as.character)

    print("Finished melting correlation matrix")

    cor_pvals <- reshape2::melt(cormat[[3]])
    cor_pvals[,1:2] <- apply(cor_pvals[,1:2],2,as.character)
    print("Finished melting p-value matrix")
#        str(cor_vals)
#        str(cor_pvals)
        cor_vals$pvals <- cor_pvals[,3]

        sub_vals <- dplyr::filter(cor_vals,pvals <=0.1)
        out_list <- list()
    print("Finished melting")
    
        for(i in genes_of_interest){
            sub <- dplyr::filter(sub_vals, Var1==i)
            sub <- sub_vals
            sub$qvals <- p.adjust(sub$pvals,method="BH")
            self <- data.frame(cbind(sub[1,1],sub[1,1],1,0,0))
            colnames(self) <- colnames(sub)
#            str(self)
#           str(sub)
            sub <- rbind(sub,self,stringsAsFactors=F)
            sub[3:5] <- apply(sub[3:5],2,as.numeric)
            out_list[[i]] <- dplyr::filter(sub, qvals <=0.1)
            print(paste0("Finished ",i))
        }

        out_list <- unique(do.call("rbind",out_list))
    str(out_list)
    colnames(out_list) <- c("Gene1","Gene2","Correlation","P_value","Q_value")
    out_list[,1] <- as.character(out_list[,1]); out_list[,2] <- as.character(out_list[,2])
    out_list$Quantile <- ecdf(abs(out_list$Correlation))(abs(out_list$Correlation))
    out <- dplyr::filter(out_list, Quantile>corr_cutoff)

    return(list(cormat,out_list,out,length(index)))}

get_druggable_genes <- function(coessential_gene_output,DGIDB_df,Achilles_df,toxic_cutoff=-0.5){
    require(dplyr)
    require(reshape2)
    DGIDB_df$gene_essentiality <- 0

    for(i in unique(DGIDB_df$gene_name)){ index <- grep(paste0("^",i,"$"), DGIDB_df[,1]); if(i %in% rownames(Achilles_df)){ DGIDB_df[index,"gene_essentiality"] <- median(as.matrix(Achilles_df)[i,],na.rm=TRUE)}
                                      }
    DGIDB_df <- DGIDB_df %>% group_by(drug_name) %>% mutate(Potential_toxic=ifelse(any(gene_essentiality <= toxic_cutoff),TRUE,FALSE)) %>% data.frame
    

    print("Getting druggable genes")
    if(is.data.frame(coessential_gene_output)==FALSE){
        df <- coessential_gene_output[[3]]} else {df <- coessential_gene_output}
##    str(df)
    df$Targetable <- df$Gene2 %in% DGIDB_df[,1]
    drug_list <- lapply(sort(unique(DGIDB_df$gene_name)), function(x) unique(toupper(dplyr::filter(DGIDB_df, gene_name==x)$drug_name)))
    tox_list <- lapply(sort(unique(DGIDB_df$gene_name)), function(x) dplyr::filter(DGIDB_df, gene_name==x)$Potential_toxic)

    print("Working on potential toxicity")


    names(drug_list) <- sort(unique(DGIDB_df$gene_name))
    names(tox_list) <- sort(unique(DGIDB_df$gene_name))
    df$Drugs <- unlist(lapply(df$Gene2, function(x) paste0(drug_list[[x]],collapse="~")))
    df$Toxic <- unlist(lapply(df$Gene2, function(x) paste0(tox_list[[x]],collapse="~")))
    df$Toxic_fraction <- unlist(lapply(df$Gene2, function(x) sum(as.numeric(tox_list[[x]]))/(length(tox_list[[x]])+1e-99)))

    druggable_list <- list(df, drug_list, tox_list)
    
    return(druggable_list)}

get_final_candidate_scores <- function(druggable_genes_output, TF_score_output, essentiality_df, CGC_list, TF_list){
    require(dplyr)

    
    #######Getting essentiality scores for Genes, obtaining p-value of median essentiality relative to others#########
    druggable_genes_output$Essentiality <- unlist(lapply(druggable_genes_output$Gene2, function(x) median(unlist(essentiality_df[x,]),na.rm=TRUE)))
    druggable_genes_output$Essentiality_pval <- 2*(1-pnorm(abs(druggable_genes_output$Essentiality),mean(druggable_genes_output$Essentiality,na.rm=TRUE),sd(druggable_genes_output$Essentiality,na.rm=TRUE)))
    druggable_genes_output$Always_essential <- ifelse(druggable_genes_output$Essentiality <0 & druggable_genes_output$Essentiality_pval <0.05,TRUE,FALSE)

########Get TF score clusters for each gene #####################
    druggable_genes_output$Cluster_Gene1 <- unlist(lapply(druggable_genes_output$Gene1,function(x) paste0(unique(TF_score_output[grep(paste0("^",x,"-|-",x,"-|-",x,"$"), TF_score_output$Cluster_drivers),"Disease"]),collapse="|")))
    druggable_genes_output$Cluster_Gene2 <- unlist(lapply(druggable_genes_output$Gene2,function(x) paste0(unique(TF_score_output[grep(paste0("^",x,"-|-",x,"-|-",x,"$"), TF_score_output$Cluster_drivers),"Disease"]),collapse="|")))

    druggable_genes_output$Cluster_TF1 <- unlist(lapply(druggable_genes_output$Cluster_Gene1, function(x) paste0(unique(dplyr::filter(TF_score_output, Disease %in% unlist(strsplit(x,"\\|")))$Cluster_drivers),collapse="|")))
    druggable_genes_output$Cluster_TF2 <- unlist(lapply(druggable_genes_output$Cluster_Gene2,function(x) paste0(unique(dplyr::filter(TF_score_output, Disease %in% unlist(strsplit(x,"\\|")))$Cluster_drivers),collapse="|")))

    
    
####### Get number of drugs gene recurrence and number of clusters per gene ########
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Num_Clusters=length(unique(unlist(strsplit(Cluster_Gene1, "~"))))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Drugs) %>% mutate(Num_Drugs=length(unique(unlist(strsplit(Drugs,"~"))))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene2) %>% mutate(Gene_recurrence=length(unique(Gene1))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Targetable=(Num_Drugs>=1)) %>% data.frame

############### Check if Gene is in CGC list or a TF ########
    druggable_genes_output$Is_TF <- druggable_genes_output$Gene2 %in% TF_list
    druggable_genes_output$CGC <- druggable_genes_output$Gene2 %in% CGC_list

########### Get Mean rank for TF in clusters where it's important #####
    druggable_genes_output$TF_score_rank <- unlist(lapply(1:nrow(druggable_genes_output), function(x) mean(unlist(lapply(strsplit(druggable_genes_output[x,"Cluster_TF1"],"\\|"), function(y) lapply(strsplit(y,"-"), function(z) grep(druggable_genes_output[x,1],z)))))))
    
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Simple_score=Targetable*((abs(Correlation))+Essentiality+1/Num_Clusters)) %>% mutate(Simple_score_with_TF_rank=Simple_score+1/TF_score_rank) %>% mutate(Simple_score_TF_rank_is_TF=Simple_score+1/TF_score_rank+Is_TF) %>% mutate(Simple_score_TF_rank_is_TF_is_CGC=Simple_score+1/TF_score_rank+Is_TF+CGC) %>% mutate(Simple_score_TF_rank_is_TF_is_CGC_Tox=Simple_score+1/TF_score_rank+Is_TF+CGC+1-Toxic_fraction) %>% data.frame

return(druggable_genes_output)

}


profile_environment <- function(obj_env=get(ls()),munit="Gb"){

    sizes <- unlist(lapply(obj_env, function(x) format(object.size(get(x)),units=munit)))
    df <- data.frame(cbind(obj_env,sizes),stringsAsFactors=F)

    colnames(df) <- c("Object","Sizes")
    df[,2] <- as.numeric(gsub(paste0(" ",munit),"",df[,2]))
    df$Unit <- munit
    df <- df[order(df$Sizes),]

    return(df)
}

funseq_to_granges <- function(funseq_input, format=c("bed","vcf"),header=TRUE,strand=NULL){
    require(data.table)
    require(GenomicRanges)

    if(format=="bed"){

        infile <- as.data.frame(fread(funseq_input,header=header, sep='\t', stringsAsFactors=F))
        if(nrow(infile)<=1){
            print("Doing Nothing")
         } else{

        colnames(infile)[1] <- gsub("#","",colnames(infile)[1])
        str(infile)
        index <- as.numeric(ncol(infile))
        
        metadata <- infile[,index]
        
        
        metadata <- lapply(metadata, function(x) strsplit(as.character(x),";"))
        
        int <- do.call("rbind",do.call("rbind", metadata))
        print(dim(int))
        str(int)
                                        #        str(int)
        if(header==TRUE){
            cols <- unlist(strsplit(colnames(infile)[index],";"))} else {
                cols <- c("gerp","cds","variant.annotation.cds","network.hub","gene.under.negative.selection","ENCODE.annotated","hot.region","motif.analysis","sensitive","ultra.sensitive","ultra.conserved","target.gene[known_cancer_gene/TF_regulating_known_cancer_gene/differential_expressed_in_cancer,actionable_gene]","user.annotations","coding.score","noncoding.score","recurrence")}
        print(length(cols))
        print(cols)
        colnames(int) <- cols[1:ncol(int)]
##        str(int)    
        infile[,index] <- NULL

        int2 <- cbind(infile,int,stringsAsFactors=F)
##        str(int2)
        if(is.null(strand)){
            strand_index <- 3

            colnames(int2)[strand_index+1:strand_index+3] <- c("Ref","Alt","Sample")
        str(int2)
       
        
        gr <- GRanges(seqnames=Rle(int2[,1]), ranges= IRanges(start=int2[,2], end=int2[,3]),strand="*")

    } else{
        print("working 2")
        strand_index <- grep("^strand$",colnames(int2),ignore.case=TRUE)

        colnames(int2)[strand_index+1:strand_index+3] <- c("Ref","Alt","Sample")
##        str(int2)
        gr <- GRanges(seqnames=Rle(int2[,1]), ranges= IRanges(start=int2[,2], end=int2[,3]),strand=int2[strand_index])

        }


        gr@elementMetadata@listData <- as.list(int2[(strand_index+1):ncol(int2)])

        gr$coding.score <- as.numeric(gr$coding.score)
        gr$noncoding.score <- as.numeric(gr$noncoding.score)




    }} else { print("Not available yet")}

                                            return(gr)
}


get_funseq_recurrence <- function(output_dir,outfile,extension=".bed") {

    require(dplyr)
    require(reshape2)

    dir_list <- list.files(output_dir ,"TCGA")
    file_list <- list.files(paste0(output_dir,dir_list),extension)
    samples <- unique(gsub("_candidates.bed","",file_list))
    str(samples)
    system(paste0("cat ", output_dir,samples[1],"/",grep(samples[1],file_list,value=TRUE), " > ",outfile))
    for(i in 2:length(samples)){
        system(paste0("tail -n +2 ",paste0(output_dir,samples[i],"/",grep(samples[i],file_list,value=TRUE)), " > int"))
        system(paste0("cat int >> ", outfile))}
    print(paste0("Finished concatenating ", length(samples)," files"))

    infile <- read.delim(outfile, header=TRUE, sep='\t', stringsAsFactors=F)
    infile2 <- infile %>% group_by(Chrom, Start, End) %>% mutate(N_with_Alteration=length(unique(Sample))) %>% mutate(Mutated_samples=paste0(Sample,collapse="|")) %>% data.frame
    infile2 <- infile2 %>% group_by(Chrom, Start, End,Alt) %>% mutate(N_with_SNV=length(unique(Sample))) %>% mutate(Samples_with_SNV=paste0(Sample, collapse="|")) %>% data.frame
    str(infile2)
    write.table(infile2, outfile, quote=FALSE, sep='\t',row.names=FALSE)
    system("rm int")
        
}

bed_to_granges_dynamic <- function(bed,header=FALSE,nlines=-1){
    require(GenomicRanges)

    
    bed <- read.table(bed, stringsAsFactors=F, sep="\t",header=header,nrows=nlines)
##    str(bed)
    if(nrow(bed) <=1 && header==TRUE){ print("Doing Nothing")} else if(nrow(bed) <1 && header==FALSE){ print("Doing nothing, empty file")} else if (nrow(bed) ==1 && header==FALSE){
    gr <- GRanges(seqnames=Rle(bed[,1]), ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))
    gr@elementMetadata@listData <- as.list(bed[4:ncol(bed)])

    return(gr)}
                                                                else{
                                                                    gr <- GRanges(seqnames=Rle(bed[,1]), ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])))
    gr@elementMetadata@listData <- as.list(bed[4:ncol(bed)])

                                                                    return(gr)}}

bedpe_to_granges <- function(bedpe,ID_col=NULL,ranges_1=1:3,ranges_2=4:6,strand=NULL,header=FALSE,filter_trans=FALSE,el_names="Link_"){

    if(is.character(bedpe)){
    bedpe <- read.delim(bedpe, stringsAsFactors=F, sep='\t',header=header)} else if(is.data.frame(bedpe)){ bedpe <- bedpe}

    if(filter_trans==TRUE){ bedpe <- bedpe[which(bedpe[,ranges_1[1]] ==bedpe[,ranges_2[1]]),]
                        print("Finished filtering trans interactions") } else{ bedpe <- bedpe}

##    str(bedpe)
    metadata_cols <- setdiff(1:ncol(bedpe),c(ranges_1,ranges_2))

    table1 <- bedpe[c(ranges_1,metadata_cols)]
    table2 <- bedpe[c(ranges_2,metadata_cols)]
   ## print(nrow(bedpe))
    link_names <- if(!is.null(ID_col)){
       link_names <- bedpe[,ID_col]} else { link_names <-  paste0(el_names,1:nrow(bedpe))}

   ## str(table1)
   ## str(table2)
   ## str(link_names)

    colnames(table1)[1:3] <- c("Chr","Start","End")
    colnames(table2)[1:3] <- c("Chr","Start","End")

    table1$Name <- link_names
    table2$Name <- link_names

    table_out <- rbind(table1, table2, stringsAsFactors=F)
##    str(table_out)

    gr_out <- table_to_granges(table_out)
    ##print(gr_out)
    out_list <- split(gr_out,gr_out$Name)
    

##    out_list <- list();
##    for(i in 1:length(link_names)){
##        if(!is.null(strand)){ gr1 <- table_to_granges(table1[i,],strand=last(ranges_1))
##                              gr2 <-table_to_granges(table2[i,],strand=last(ranges_1))} else{
##                                  gr1 <- table_to_granges(table1[i,])
##                                  gr2 <- table_to_granges(table2[i,])}

                                    ##print(gr1)
                                    ##print(gr2)
##                                    gr <- union(gr1,gr2)
##                                    gr$Name <- link_names[i]
                                    ##print(gr)
##                                    out_list[[link_names[[i]]]] <- gr
##                                }

##    out_list <- GRangesList(out_list)
    return(out_list)
    
}


table_to_granges_bedpe <- function(table,ID_col=NULL,ranges_1=1:3,ranges_2=4:6,strand=NULL,header=FALSE,filter_trans=FALSE){
    bedpe <- table

    if(filter_trans==TRUE){ bedpe <- bedpe[which(bedpe[,ranges_1[1]] ==bedpe[,ranges_2[1]]),]
                        print("Finished filtering trans interactions") } else{ bedpe <- bedpe}

##    str(bedpe)
    metadata_cols <- setdiff(1:ncol(bedpe),c(ranges_1,ranges_2))

    table1 <- bedpe[c(ranges_1,metadata_cols)]
    table2 <- bedpe[c(ranges_2,metadata_cols)]
   ## print(nrow(bedpe))
    link_names <- if(!is.null(ID_col)){
       link_names <- bedpe[,ID_col]} else { link_names <-  paste0("Link_",1:nrow(bedpe))}

   ## str(table1)
   ## str(table2)
   ## str(link_names)

    colnames(table1)[1:3] <- c("Chr","Start","End")
    colnames(table2)[1:3] <- c("Chr","Start","End")

    table1$Name <- link_names
    table2$Name <- link_names

    table_out <- rbind(table1, table2, stringsAsFactors=F)
##    str(table_out)

    gr_out <- table_to_granges(table_out)
    ##print(gr_out)
    out_list <- split(gr_out,gr_out$Name)
    

##    out_list <- list();
##    for(i in 1:length(link_names)){
##        if(!is.null(strand)){ gr1 <- table_to_granges(table1[i,],strand=last(ranges_1))
##                              gr2 <-table_to_granges(table2[i,],strand=last(ranges_1))} else{
##                                  gr1 <- table_to_granges(table1[i,])
##                                  gr2 <- table_to_granges(table2[i,])}

                                    ##print(gr1)
                                    ##print(gr2)
##                                    gr <- union(gr1,gr2)
##                                    gr$Name <- link_names[i]
                                    ##print(gr)
##                                    out_list[[link_names[[i]]]] <- gr
##                                }

##    out_list <- GRangesList(out_list)
    return(out_list)
    
}


sort_gr <- function(gr){
    require(GenomicRanges)
    gr <- sort(sortSeqlevels(gr))

    return(gr)}

get_unique_ranges <- function(gr, column_name){
    grl <- lapply(unique(elementMetadata(gr)[,column_name]), function(x) gr[which(elementMetadata(gr)[,column_name] ==x)])

    names(grl) <- unique(elementMetadata(gr)[,column_name])
                                        #    print(grl)
    grl_int <- lapply(grl, reduce)
    names(grl_int) <- names(grl)

#    print(grl[unique(elementMetadata(gr)[,column_name])][[1]])
#    print(grl[c(1,2)])
    grl2 <- lapply(names(grl_int), function(x) setdiff_with_metadata(grl[[x]], unlist(GRangesList(grl[setdiff(names(grl),x)]))))
    
    names(grl2) <- names(grl)
    grl_final <- GRangesList(lapply(names(grl2), function(x) intersect_with_metadata(grl[[x]], grl2[[x]])))
    names(grl_final) <- names(grl)
    return(grl_final)
}

motif_analysis_to_matrix <- function(funseq_gr){
    require(stringi)
    ##Split motif column into separate columns ##
    
    funseq_gr$Motif_Gain <- unlist(lapply(funseq_gr$motif.analysis, function(x)  unlist(strsplit(x,"MOTIFG="))[2]))


    funseq_gr$Motif_Break <- gsub("MOTIFBR=","",unlist(lapply(funseq_gr$motif.analysis, function(x)  unlist(strsplit(x,"MOTIFG="))[1])))

    #print(funseq_gr)

    print("Finished Parsing Motif analysis")


    ##Initialize for loop ##
    
    gr_ranges <- reduce(funseq_gr)
    print(length(unique(funseq_gr$Sample)))
    out_mat <- matrix(0, nrow=length(unique(funseq_gr$Sample)))

    rownames(out_mat) <- unique(funseq_gr$Sample)
    out_mat_gained <- out_mat

    total_broken <- c()
    total_gained <- c()
    print("Making matrix")
                      
    for(i in 1:length(gr_ranges)){
        sub_gr <- reduce(gr_ranges[i])
        df <- unique(as.data.frame(elementMetadata(intersect_with_metadata(funseq_gr, sub_gr))))
        


       for(k in unique(df$motif.analysis)){
            sub_df <- df[which(df$motif.analysis ==k),]
            samples <- unique(sub_df$Sample)

            motifs_broken <- unique(unlist(strsplit(unique(sub_df$Motif_Break),",")))
            motifs_gained <- unique(unlist(strsplit(unique(sub_df$Motif_Gain),",")))

            
            motifs_broken <- gsub("\\.","",unique(unlist(lapply(motifs_broken, function(x) unlist(strsplit(x,"#"))[1]))))
            motifs_gained <- gsub("\\.","",unique(unlist(lapply(motifs_gained, function(x) unlist(strsplit(x,"#|_"))[1]))))

            motifs_broken <- motifs_broken[motifs_broken != ""]
            motifs_gained[is.na(motifs_gained)] <- ""
            motifs_gained <- motifs_gained[motifs_gained != ""]
            total_broken <- unique(c(total_broken,motifs_broken))
            total_gained <- unique(c(total_gained,motifs_gained))

            
            ##Working Motifs Broken ###
            if(length(motifs_broken) >=1){
            sub_mat <- matrix(0,nrow=nrow(out_mat),ncol=length(motifs_broken))
            dimnames(sub_mat) <-  list(rownames(out_mat),motifs_broken)
            sub_mat[samples,motifs_broken] <- 1
        
            
            mcol <- match(colnames(sub_mat),colnames(out_mat))
            if(all(is.na(mcol)==TRUE)){
#                print("Working 1")
                out_mat <- merge(out_mat, sub_mat,by=0,all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]



#            } else if (all(is.na(mcol) ==FALSE)) {


                
#                print("Working 2")
#                out_mat[,mcol] <- out_mat[,mcol]+sub_mat
               
            } else {
#                print("Working 3")
                out_mat[,intersect(colnames(out_mat),colnames(sub_mat))] <- out_mat[,intersect(colnames(out_mat),colnames(sub_mat))]+sub_mat[,intersect(colnames(out_mat),colnames(sub_mat))]
                out_mat <- merge(out_mat, sub_mat[,setdiff(colnames(sub_mat),colnames(out_mat))],by=0, all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                colnames(out_mat)[grep("y", colnames(out_mat))] <- setdiff(colnames(sub_mat),colnames(out_mat))
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]


            }}

                
#########################################################################################          
            ##Working Motifs Gained ##
        
                 if(length(motifs_gained) >=1){
            sub_mat_gained <- matrix(0,nrow=nrow(out_mat_gained),ncol=length(motifs_gained))
            dimnames(sub_mat_gained) <-  list(rownames(out_mat_gained),motifs_gained)
            sub_mat_gained[samples,motifs_gained] <- 1


        

            mcol <- match(colnames(sub_mat_gained),colnames(out_mat_gained))
            if(all(is.na(mcol)==TRUE)){
                out_mat_gained <- merge(out_mat_gained, sub_mat_gained,by=0,all=TRUE)
                names_of_rows  <- out_mat_gained[,"Row.names"]
                names_of_cols <- colnames(out_mat_gained)
                out_mat_gained$Row.names <- NULL
                out_mat_gained <- as.matrix(out_mat_gained)
                out_mat_gained <- apply(out_mat_gained,2,as.numeric)
                rownames(out_mat_gained) <- names_of_rows
                colnames(out_mat_gained) <- names_of_cols[2:length(names_of_cols)]




#            } else if (all(is.na(mcol)==FALSE)){
#
#                out_mat_gained[,mcol] <- out_mat_gained[,mcol]+sub_mat_gained

            
            } else {
                                        #                print("Working 3")
                 out_mat_gained[,intersect(colnames(out_mat_gained),colnames(sub_mat_gained))] <- out_mat_gained[,intersect(colnames(out_mat_gained),colnames(sub_mat_gained))]+sub_mat_gained[,intersect(colnames(out_mat_gained),colnames(sub_mat_gained))]
                out_mat_gained <- merge(out_mat_gained, sub_mat_gained[,setdiff(colnames(sub_mat_gained),colnames(out_mat_gained))],by=0, all=TRUE)
                names_of_rows  <- out_mat_gained[,"Row.names"]
                colnames(out_mat_gained)[grep("y", colnames(out_mat_gained))] <- setdiff(colnames(sub_mat_gained),colnames(out_mat_gained))
                names_of_rows  <- out_mat_gained[,"Row.names"]
                names_of_cols <- colnames(out_mat_gained)
                out_mat_gained$Row.names <- NULL
                out_mat_gained <- as.matrix(out_mat_gained)
                out_mat_gained <- apply(out_mat_gained,2,as.numeric)
                rownames(out_mat_gained) <- names_of_rows
                colnames(out_mat_gained) <- names_of_cols[2:length(names_of_cols)]


             }} else { print("No motifs altered")}
        
        }}

    
   combined_mcol <- match(colnames(out_mat),colnames(out_mat_gained))
  
    output <- list(out_mat[,2:ncol(out_mat)], out_mat_gained[,2:ncol(out_mat_gained)], add_non_conformable_matrices(out_mat[,2:ncol(out_mat)], out_mat_gained[,2:ncol(out_mat_gained)]))
    names(output) <- c("Motifs_Lost","Motifs_Gained","Sum_Affected")
                   
    return(output)

}


add_non_conformable_matrices <- function(mat1, mat2){

    

    union_cols <- union(colnames(mat1),colnames(mat2))
    union_rows <- union(rownames(mat1), rownames(mat2))
    union_mat <- matrix(0,nrow=length(union_rows),ncol=length(union_cols))
    dimnames(union_mat) <- list(union_rows, union_cols)

    union_mat[rownames(mat1),colnames(mat1)] <- mat1[rownames(mat1),colnames(mat1)]
    union_mat[rownames(mat2),colnames(mat2)] <- union_mat[rownames(mat2),colnames(mat2)]+mat2[rownames(mat2),colnames(mat2)]
    return(union_mat)
    }
        
network_analysis_to_matrix <- function(funseq_gr){
    gr_ranges <- reduce(funseq_gr)
    print(length(unique(funseq_gr$Sample)))
    out_mat <- matrix(0, nrow=length(unique(funseq_gr$Sample)))

    rownames(out_mat) <- unique(funseq_gr$Sample)
    for(i in 1:length(gr_ranges)){
        sub_gr <- reduce(gr_ranges[i])
        df <- unique(as.data.frame(elementMetadata(intersect_with_metadata(funseq_gr, sub_gr))))

        for(k in unique(df$network.hub)){
            sub_df <- df[which(df$network.hub ==k),]
            samples <- unique(sub_df$Sample)
            genes_affected <- unique(unlist(strsplit(unique(sub_df$network.hub),",")))
            
            genes_affected <- gsub("\\.","",unique(unlist(lapply(genes_affected, function(x) unlist(strsplit(x,":"))[1]))))
            genes_affected <- genes_affected[genes_affected != ""]
            ##Working Genes Affected ###
            if(length(genes_affected) >=1){
            sub_mat <- matrix(0,nrow=nrow(out_mat),ncol=length(genes_affected))
            dimnames(sub_mat) <-  list(rownames(out_mat),genes_affected)
            sub_mat[samples,genes_affected] <- 1


            mcol <- match(colnames(sub_mat),colnames(out_mat))
            if(all(is.na(mcol)==TRUE)){
#                print("Working 1")
                out_mat <- merge(out_mat, sub_mat,by=0,all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]



#            } else if (all(is.na(mcol) ==FALSE)) {


                
#                print("Working 2")
#                out_mat[,mcol] <- out_mat[,mcol]+sub_mat
               
            } else {
#                print("Working 3")
                out_mat[,intersect(colnames(out_mat),colnames(sub_mat))] <- out_mat[,intersect(colnames(out_mat),colnames(sub_mat))]+sub_mat[,intersect(colnames(out_mat),colnames(sub_mat))]
                out_mat <- merge(out_mat, sub_mat[,setdiff(colnames(sub_mat),colnames(out_mat))],by=0, all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                colnames(out_mat)[grep("y", colnames(out_mat))] <- setdiff(colnames(sub_mat),colnames(out_mat))
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]


            }} else { print("No motifs altered")}}
    }

    return(out_mat[,2:ncol(out_mat)])
}
 
                                             
gene_analysis_to_matrix <- function(funseq_gr){
    gr_ranges <- reduce(funseq_gr)
    print(length(unique(funseq_gr$Sample)))
    out_mat <- matrix(0, nrow=length(unique(funseq_gr$Sample)))

    rownames(out_mat) <- unique(funseq_gr$Sample)
    for(i in 1:length(gr_ranges)){
        sub_gr <- reduce(gr_ranges[i])
        df <- unique(as.data.frame(elementMetadata(intersect_with_metadata(funseq_gr, sub_gr))))

        for(k in unique(df$target.gene.known_cancer_gene.TF_regulating_known_cancer_gene.differential_expressed_in_cancer.actionable_gene.)){

            sub_df <- df[which(df$target.gene.known_cancer_gene.TF_regulating_known_cancer_gene.differential_expressed_in_cancer.actionable_gene. ==k),]

            samples <- unique(sub_df$Sample)
            
            genes_affected <- gsub("\\s*\\([^\\)]+\\)","",unique(sub_df$target.gene.known_cancer_gene.TF_regulating_known_cancer_gene.differential_expressed_in_cancer.actionable_gene.))
            genes_affected <- unique(unlist(strsplit(unique(genes_affected),",")))
            genes_affected <- gsub("\\s*\\[[^\\]]+\\]","",genes_affected)
            genes_affected <- unlist(lapply(genes_affected,function(x) unlist(strsplit(x,"\\["))[1]))
            genes_affected <- unlist(lapply(genes_affected,function(x) gsub("\\[|\\]","",x)))            

            genes_affected <- genes_affected[genes_affected != ""]

            ##Working Genes Affected ###
            if(length(genes_affected) >=1){
            sub_mat <- matrix(0,nrow=nrow(out_mat),ncol=length(genes_affected))
            dimnames(sub_mat) <-  list(rownames(out_mat),genes_affected)
            sub_mat[samples,genes_affected] <- 1


            mcol <- match(colnames(sub_mat),colnames(out_mat))
            if(all(is.na(mcol)==TRUE)){
#                print("Working 1")
                out_mat <- merge(out_mat, sub_mat,by=0,all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]



#            } else if (all(is.na(mcol) ==FALSE)) {


                
#                print("Working 2")
#                out_mat[,mcol] <- out_mat[,mcol]+sub_mat
               
            } else {
#                print("Working 3")
                out_mat[,intersect(colnames(out_mat),colnames(sub_mat))] <- out_mat[,intersect(colnames(out_mat),colnames(sub_mat))]+sub_mat[,intersect(colnames(out_mat),colnames(sub_mat))]
                out_mat <- merge(out_mat, sub_mat[,setdiff(colnames(sub_mat),colnames(out_mat))],by=0, all=TRUE)
                names_of_rows  <- out_mat[,"Row.names"]
                colnames(out_mat)[grep("y", colnames(out_mat))] <- setdiff(colnames(sub_mat),colnames(out_mat))
                names_of_cols <- colnames(out_mat)
                out_mat$Row.names <- NULL
                out_mat <- as.matrix(out_mat)
                out_mat <- apply(out_mat,2,as.numeric)
                rownames(out_mat) <- names_of_rows
                colnames(out_mat) <- names_of_cols[2:length(names_of_cols)]


            }} else { print("No genes altered")}}
    }

    return(out_mat[,2:ncol(out_mat)])
}


binarize_matrix <- function(mat, cutoff){
    names_of_rows <- rownames(mat)
    mat <- apply(mat,2, function(x) as.numeric(x>=cutoff))
    rownames(mat) <- names_of_rows
    return(mat)}

get_true_peak_gene_links <- function(P_G_link_mat,peak_height_mat,cutoff){
    P_G_link_mat <- P_G_link_mat[intersect(colnames(P_G_link_mat), colnames(peak_height_mat))]
    peak_height_mat <- peak_height_mat[intersect(colnames(P_G_link_mat), colnames(peak_height_mat))]
        
    peaks <- unlist(lapply(rownames(P_G_link_mat), function(x) unlist(strsplit(x,"~"))[1]))
##    str(peaks)
    sub_peak_mat <- peak_height_mat[peaks,]
##    str(sub_peak_mat)
    binary_sub <- binarize_matrix(sub_peak_mat, cutoff)

    out <- P_G_link_mat*binary_sub
    return(out)
}
summarize_binary_matrix <- function(mat,metadata,group_column,sample_column){
    require(dplyr)
    require(reshape2)
    df <- reshape::melt(mat)
    df[,1:2] <- apply(df[,1:2],2, as.character)
    colnames(df) <- c("Sample","Gene","value")
    df[,"Metadata"] <- metadata[df$Sample,group_column]
    df2 <- df %>% group_by(Gene,Metadata) %>% summarise(Count=sum(value)) %>% data.frame

    return(df2)}


compare_clustering <- function(df1, df2,TF_column=c("Cluster_drivers","Signif_drivers"),df2_description=NULL){
    require(dplyr)
    df2 <- df2[df1$Sample,]
    novel <- lapply(df1$Sample, function(x) paste0(setdiff(unlist(strsplit(df1[x,TF_column],"-")),unlist(strsplit(df2[x,TF_column],"-"))),collapse="-"))
    
    common <- lapply(df1$Sample, function(x) paste0(intersect(unlist(strsplit(df1[x,TF_column],"-")),unlist(strsplit(df2[x,TF_column],"-"))),collapse="-"))

    names(novel) <- df1$Sample
    names(common) <- df1$Sample    

    
    df_sub <- merge(df1,df2[,unique(c(colnames(df2), "Disease"))],by="Sample")
    df_sub$Novel_TFs <- unlist(novel[df_sub$Sample])
    df_sub$Common_TFs <- unlist(common[df_sub$Sample])
    
    col_index <- grep("\\.y", colnames(df_sub))
    colnames(df_sub)[col_index] <- paste0(df2_description,"_", colnames(df_sub)[col_index])
    colnames(df_sub) <- gsub("\\.y|\\.x", "", colnames(df_sub))

    df_sub <- df_sub %>% group_by(Sample) %>% mutate(n_novel=length(unlist(strsplit(Novel_TFs,"-")))) %>% mutate(n_common=length(unlist(strsplit(Common_TFs,"-"))))  %>% data.frame

    
    return(df_sub)
}

zscore_df <- function(df, margin=c("row","column"),central_tendency=c("mean","median","median_deviation")){
    if(length(central_tendency)>1){ central_tendency <- "mean"}
    if(length(margin)>1){ margin <- "row"}
    if(margin == "row"){
        if(central_tendency=="mean"){df <- t(apply(df,1,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))} else if(central_tendency=="median"){ df <- t(apply(df,1,function(x) (x-median(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))} else if (central_tendency=="median_deviation"){ df <- t(apply(df,1,function(x) (x-median(x,na.rm=TRUE))))}

    } else if(margin=="column") { if(central_tendency=="mean"){df <- apply(df,2,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))} else if (central_tendency=="median"){ df <- apply(df,2,function(x) (x-median(x,na.rm=TRUE))/sd(x,na.rm=TRUE))} else if (central_tendency=="median_deviation"){ df <- apply(df,2,function(x) (x-median(x,na.rm=TRUE)))}
                                                                                                                                                                                                                                                         }
    df[is.na(df)] <- 0
    df <- as.data.frame(df)

    return(df)}

set_highlight <- function(df,column_to_search,search_value,column_to_change,new_value, non_highlight_val=NULL){
    
    if(length(search_value) > 1){
        highlight_index <- which(df[,column_to_search] %in% search_value)} else {
            highlight_index <- which(df[,column_to_search] == search_value) }

    df[highlight_index,column_to_change] <- new_value
    non_highlight <- setdiff(1:nrow(df), highlight_index)
    
    if(is.null(non_highlight_val) == TRUE){
        df[non_highlight,column_to_change] <- "#D3D3D303"} else{ df[non_highlight,column_to_change] <- non_highlight_val}

    return(df)}
drug_penalty <- function(vector) {
    index <- which(vector ==0)
    penalty <- unlist(lapply(vector, function(x) abs((2/x^2)*exp(1+x/max(x))-x)))
    penalty[index] <- 10000
    return(penalty)}
    
liftover_dir <- function(dir,extension,genome,header=FALSE, liftover_file,filetype=c("vcf","bed","bedpe")) {
    require(VariantAnnotation)
    require(rtracklayer)

    file_list <- list.files(dir,extension)
    for(i in file_list){
        if(filetype=="vcf"){
        vcf <- readVcf(paste0(dir,i),genome)
        gr <- rowRanges(vcf)
        gr$paramRangeID <- NULL
        gr$FILTER <- NULL
        gr$QUAL <- NULL

        sample_id <- toupper(unlist(strsplit(i,"--"))[1])
        gr$Sample <- sample_id} else if(filetype=="bed"){
            gr <- bed_to_granges_dynamic(paste0(dir,i),header=header)} else if(filetype=="bedpe"){
            gr <- unlist(bedpe_to_granges(paste0(dir,i),NULL,header=header))}
        gr_lift <- unlist(liftOver(gr,import.chain(liftover_file)))
        gr_to_bed(gr_lift,paste0(dir,gsub(extension,paste0("_",genome,".bed"),i)), TRUE, FALSE,FALSE,TRUE)
        print(paste0("Finished ", grep(i, file_list)," of ", length(file_list)))
    }}
    
scale_TF_score <- function(color_list=NULL,extension=".txt",scales=c("fill","color"),...){
    require(ggplot2)
##    print(color_list)
    if(is.null(color_list)==TRUE){
##        print("1")
        cols <- readRDS("~/color_scale_TF_score_vec.rds")} else if(is.character(color_list) ==TRUE & !is.vector(color_list)) { if(extension == ".rds") { print("reading RDS"); cols <- readRDS(color_list)} else if(extension==".txt"){ print("reading txt");cols <- read.table(color_list,sep='\t',stringsAsFactors=F,header=T,comment.char="$")}
                                                                                                      if(is.data.frame(cols) ==TRUE){ df <- cols; cols <- df[,2]; names(cols) <- df[,1]} else{ cols <- cols}} else if(is.vector(color_list)==TRUE){ cols <- color_list}
                                       
    
    g <- ggplot2:::manual_scale(scales, values=cols)
    return(g)
}
run_pca <- function(df, num_components=20, center=TRUE,scale=TRUE){
    require(irlba)
    comp_run <- prcomp_irlba(df,n=num_components,center=center,scale=scale)
    components <- as.data.frame(comp_run$rotation)
#    str(components)
    rownames(components) <- colnames(df)
    components$Sample <- colnames(df)
    colnames(components) <- gsub("PC","Dim",colnames(components))
    importances <- as.data.frame(t(summary(comp_run)$importance))
    return(list(components,importances))
}
umap_tuning <- function(df,n_neighbors_vec,spread_vec,seed=20,minPoints=4) {
    out_df <- data.frame(stringsAsFactors=F)
    for(i in n_neighbors_vec){
        print(paste0("Working ",i))
        sub_list <- lapply(spread_vec, function(x) run_umap(df,n_neighbors=i,spread=x,seed=seed))
                           for(j in 1:length(sub_list)) {
                               sub <- sub_list[[j]]
                               sub$Neighbors <- i
                               sub$Spread <- spread_vec[j]
                               print(paste0("Finished spread ",spread_vec[j]))
                               sub_list[[j]] <- sub
                           }
                           sub_list <- lapply(sub_list, function(x) get_clusters_tsne(x,method="dbscan",minPoints=minPoints,kmin_samples=3,keep_outliers=FALSE, outlier_method="Centroid"))
                           int_df <- do.call("rbind", sub_list)
                           out_df <- rbind(out_df,int_df, stringsAsFactors=F)
                           print(paste0("Finished ",i))}
                       
    return(out_df)
                       
}
k_distance_curve_optimization <- function(df, minPts_vec,kmin_points){
    require(splines)

    out_df <- data.frame(stringsAsFactors=F)
    for(i in minPts_vec){
        print(paste0("Working ",i))
        eps <- get_optimal_dbscan_eps(df,min_samples=i,kmin_samples=kmin_points,opt_count=1)
                          #print(eps)
                          knn_df <- readRDS("knn_df.rds")
                          knn_df$Distance <- round(knn_df$Distance,3)
                          knn_df$Smoothed_Distance <- smooth.spline(knn_df$Sample, knn_df$Distance,spar=0.7)$y
                          knn_df$Smoothed_Diff <- smooth.spline(knn_df$Sample, knn_df$Diff,spar=0.7)$y
        knn_df$MinPts <- i
        knn_df$KNN <- kmin_points
        knn_df$eps <- eps
        knn_df$bisection_eps <- knn_df[which(knn_df$Sample == unique(knn_df$Elbow_bisection)),"Distance"]
                          knn_df$Elbow <- round(mean(knn_df[which(knn_df$Distance==eps),"Sample"],na.rm=TRUE))
                          
                          #str(knn_df)

        out_df <- rbind(out_df,knn_df)
        print(paste0("Finished ",i))
    }
    
    return(out_df)
}

derivative <- function(x,y){
    diff <- 0
    diff2 <- 0
    for(i in 2:length(x)){
        xdiff <- x[i]-x[i-1]
        ydiff <- y[i]-y[i-1]
        diff <- c(diff, (ydiff/xdiff))

        xdiff2 <- x[i+1]-x[i]
        ydiff2 <- y[i+1]-y[i]
        diff2 <- c(diff2, (ydiff2/xdiff2))}

    diff_df <- data.frame(diff,diff2)
    diff_df$Final <- apply(diff_df,1,function(x) mean(x, na.rm=TRUE))        
#    str(diff_df)
    Final_diff <- diff_df$Final
    Final_diff[which(is.infinite(Final_diff) == TRUE)] <- 0
    Final_diff[which(is.na(Final_diff) == TRUE)] <- 0
    diff <- round(Final_diff,3)
    return(Final_diff)}
    
dist2d <- function(a,b,c) {
    ###a represents point being examined
###b represents the coordinates of the first point on line
    ###c represents the coordinates of the last point on line
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    d <- abs(det(m))/sqrt(sum(v1*v1))
    return(d)
}


get_closest <- function(coordinate1, coordinate_df){
    distances <- unlist(lapply(1:nrow(coordinate_df), function(x) dist(rbind(coordinate1, coordinate_df[x,]))))

#    str(distances)
    min_dist <- which(distances==min(distances))
    return(min_dist)}
get_furthest <- function(coordinate1, coordinate_df){
    distances <- unlist(lapply(1:nrow(coordinate_df), function(x) dist(rbind(coordinate1, coordinate_df[x,]))))

#    str(distances)
    max_dist <- which(distances==max(distances))
    return(max_dist)}

clustering_tuning <- function(tsne_df,minPoints_vec=2:10,num_distances=10, dist_interval=0.1,kmin=3,cluster_method,outliers=TRUE,outlier_method="Centroid"){
    out_df <- data.frame(stringsAsFactors=F)
    optimal_dist <- get_optimal_dbscan_eps(tsne_df,3,kmin_samples=kmin)
    distances <- sort(unique(c(seq(optimal_dist,by=-dist_interval, length.out=num_distances/2),seq(optimal_dist,by=dist_interval, length.out=num_distances/2))))
    for(i in minPoints_vec){
        int_df <- data.frame(stringsAsFactors=F)
        for(k in distances){
            tsne_df <- get_clusters_tsne(tsne_df,method="dbscan",minPoints=i,eps=k,keep_outliers=outliers,outlier_method=outlier_method)
            tsne_df$Eps <- k
            int_df <- rbind(int_df,tsne_df)}
        int_df$MinPoints <- i
        out_df <- rbind(out_df,int_df)}
    out_df$Optimal_eps <- optimal_dist
    out_df$KNN <- kmin
    return(out_df)}
    
nearest_template_prediction <- function(Achilles_df, TF_score_df,TF_score_mat,TF_column=c("Cluster_drivers","Signif_drivers","Top_5"),use_raw_score=TRUE,verbose=TRUE) {
    require(lsa)
    require(dplyr)

    predicted_template_df <- data.frame(stringsAsFactors=F)
    
    template_df <- unique(TF_score_df[c("Cluster","Disease",TF_column)])
##    gene_set <- unique(unlist(lapply(unique(TF_score_df[,TF_column]), function(x) unlist(strsplit(x,"-")))))
    gene_set <- rownames(Achilles_df)

    
    templates <- lapply(sort(unique(template_df$Disease)), function(x) unlist(strsplit(dplyr::filter(template_df, Disease==x)[,TF_column], "-")))
#    str(templates)
    templates <- lapply(templates, function(x) intersect(x, rownames(Achilles_df)))
    names(templates) <- sort(unique(template_df$Disease))

###########
    samples <- lapply(sort(unique(template_df$Disease)), function(x) dplyr::filter(TF_score_df, Disease==x)$Sample)
    names(samples) <- sort(unique(template_df$Disease))

#    str(samples)
#    str(templates)
###################################
    
    if(use_raw_score==TRUE){
        sample_template_list <- lapply(sort(unique(template_df$Disease)), function(x) rowMeans(TF_score_mat[templates[[x]],samples[[x]]],na.rm=TRUE))
        names(sample_template_list) <- sort(unique(template_df$Disease))
    }
    
    else if (use_raw_score==FALSE) {
        sample_template_list <- lapply(sort(unique(template_df$Disease)), function(x) rowMeans(TF_score_mat[templates[[x]],samples[[x]]],na.rm=TRUE)/rowMeans(TF_score_mat[templates[[x]],setdiff(colnames(TF_score_mat),samples[[x]])],na.rm=TRUE))
        names(sample_template_list) <- sort(unique(template_df$Disease))

        print("Not finished with this yet")}
#######################
#######################
    
    for(i in colnames(Achilles_df)){
        distance_list <- unlist(lapply(sort(unique(template_df$Disease)), function(x) 1-cosine(Achilles_df[templates[[x]],i],sample_template_list[[x]])))
        names(distance_list) <- sort(unique(template_df$Disease))

        random_list <- list()
    for(k in names(templates)){
        random_list[[k]] <- unlist(lapply(1:1000, function(x) 1-cosine(Achilles_df[sample(gene_set, length(templates[[k]])),i],sample_template_list[[k]])))
    }
        #lapply(random_list,function(x) print(summary(x)))

        pval_list <- unlist(lapply(names(distance_list), function(x) 1-pnorm(distance_list[[x]],mean(random_list[[x]]),sd(random_list[[x]]),lower.tail=F)))
        names(pval_list) <- names(distance_list)
        ###print(distance_list)
        ###print("Pvalue=")
        ###print(pval_list)
        candidate_template_df <- data.frame(pval_list,distance_list,names(distance_list),i,stringsAsFactors=F)
        colnames(candidate_template_df) <- c("P_val","Distance","Template","Cell_Line")
        candidate_template_df$Signif <- candidate_template_df$P_val <= 0.25
        sub_df <- dplyr::filter(candidate_template_df, Signif== TRUE)
        if(nrow(sub_df) >=1){
            candidate_template_df$Prediction <- sub_df[which(sub_df$Distance == min(sub_df$Distance)),"Template"]} else{ candidate_template_df$Prediction <- candidate_template_df[which(candidate_template_df$Distance == min(candidate_template_df$Distance)),"Template"]}

        
        predicted_template_df <- rbind(predicted_template_df,candidate_template_df)                                      
                                        #str(sub_df)
        if(verbose==TRUE){
        print(paste0("Finished ", grep(paste0("^",i,"$"), colnames(Achilles_df)), " of ",ncol(Achilles_df)))} 
    }
        return(predicted_template_df)
}


get_template_genes <- function(expression_df, tsne_df, cluster_column="Cluster",is_ranked=TRUE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,value_returned=c("Templates","Matrix")){
    require(randomForest)
    source("~/Andre_F_functions.R")
    require(dplyr)

    sub_df <- expression_df[,intersect(tsne_df$Sample, colnames(expression_df))]
    sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]

    if(method == "simple"){
        df <- zscore_df(do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample]))),"row")
        df_abs <- do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample])))
        colnames(df) <- unique(tsne_df[,cluster_column])
        colnames(df_abs) <- unique(tsne_df[,cluster_column])
        if(is_ranked == TRUE){
            
            templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=FALSE))[1:1000])
            templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:1000])
            names(templates_abs) <- colnames(df)
            names(templates_rel) <- colnames(df)
            templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            
        } else if (is_ranked==FALSE){ templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=TRUE))[1:1000])
                                                                                                                                                 templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:1000])
                                                                                                                                                 names(templates_abs) <- colnames(df)
                                                                                                                                                 names(templates_rel) <- colnames(df)

                                                                                                                                                 templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
                                                                                                                                             }
        names(templates) <- colnames(df)
        
    }

    else if(method == "randomforest"){
    
    transposed_df <- as.data.frame(t(sub_df))
    transposed_df[,cluster_column] <- as.factor(TF_score_final_TSNE[colnames(expression_df_all),cluster_column])
    uniq_classes <- unique(tsne_df[,cluster_column])
##    print(index)
##   print(cluster_column)


##    str(transposed_df)
    print("Starting Random Forest")
    rf <- randomForest(as.formula(paste0(cluster_column,"~.")),data=transposed_df,importance=TRUE)
    importance_df <- as.data.frame(importance(rf))[setdiff(colnames(importance(rf)),uniq_classes)]

    importance_df <- importance_df[order(-importance_df[,2]),]
    final_importance_matrix <- data.frame(importance(rf)[rownames(importance_df),uniq_classes],stringsAsFactors=F)
    
    colnames(final_importance_matrix) <- uniq_classes
    #str(final_importance_matrix)
    templates_abs <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, sort_order=TRUE))[1:1000])
#    str(templates_abs)
    templates_rel <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(zscore_df(final_importance_matrix),x,sort_order=TRUE))[1:1000])
#    str(templates_rel)
        ##str(final_importance_matrix)
    names(templates_abs) <- colnames(final_importance_matrix)
    names(templates_rel) <- colnames(final_importance_matrix)
    templates <- lapply(names(templates_abs), function(x) intersect(templates_abs[[x]], templates_rel[[x]])[1:template_length])


#    templates <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, decreasing=TRUE))[1:template_length])
    names(templates) <- names(templates_abs)
}
#    str(templates)

    if(value_returned=="Templates"){ return(templates)} else if(value_returned=="Matrix"){ return(final_importance_matrix)}
    
}

    

    
nearest_template_prediction_v2 <- function(reference_df,alt_df, tsne_df, cluster_column="Cluster",is_ranked=TRUE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,verbose=TRUE,random_sample=100) {
    require(lsa)
    require(dplyr)
    require(randomForest)

    predicted_template_df <- data.frame(stringsAsFactors=F)
    
    template_df <- unique(tsne_df[unique(c("Cluster",cluster_column))])
##    gene_set <- unique(unlist(lapply(unique(tsne_df[,TF_column]), function(x) unlist(strsplit(x,"-")))))
    gene_set <- rownames(reference_df)

    
    templates <- get_template_genes(reference_df, tsne_df,cluster_column,quantile_cutoff=quantile_cutoff,method=method, template_length=template_length, is_ranked=is_ranked,value_returned="Templates")
    #str(templates)
##    templates <- lapply(templates, function(x) intersect(x, rownames(reference_df)))
##    names(templates) <- sort(unique(template_df$Disease))

###########
    samples <- lapply(sort(unique(template_df[,cluster_column])), function(x) dplyr::filter(tsne_df, get(cluster_column)==x)$Sample)
    names(samples) <- sort(unique(template_df[,cluster_column]))

##    str(samples)
##    str(templates)
###################################
    
##    if(use_raw_score==TRUE){
        sample_template_list <- lapply(sort(unique(template_df[,cluster_column])), function(x) rowMeans(reference_df[templates[[x]],samples[[x]]],na.rm=TRUE))
        names(sample_template_list) <- sort(unique(template_df[,cluster_column]))
##  }
    
##    else if (use_raw_score==FALSE) {
##        sample_template_list <- lapply(sort(unique(template_df$Disease)), function(x) rowMeans(TF_score_mat[templates[[x]],samples[[x]]],na.rm=TRUE)/rowMeans(TF_score_mat[templates[[x]],setdiff(colnames(TF_score_mat),samples[[x]])],na.rm=TRUE))
##        names(sample_template_list) <- sort(unique(template_df$Disease))

##        print("Not finished with this yet")}
#######################
#######################
    
    for(i in colnames(alt_df)){
        distance_list <- unlist(lapply(sort(unique(template_df[,cluster_column])), function(x) 1-cosine(alt_df[templates[[x]],i],sample_template_list[[x]])))
        names(distance_list) <- sort(unique(template_df[,cluster_column]))

        random_list <- list()
    for(k in names(templates)){
        random_list[[k]] <- unlist(lapply(1:random_sample, function(x) 1-cosine(alt_df[sample(gene_set, length(templates[[k]])),i],sample_template_list[[k]])))

        }
        #lapply(random_list,function(x) print(summary(x)))

        pval_list <- unlist(lapply(names(distance_list), function(x) 1-pnorm(distance_list[[x]],mean(random_list[[x]],na.rm=TRUE),sd(random_list[[x]],na.rm=TRUE),lower.tail=F)))
        names(pval_list) <- names(distance_list)
        ###print(distance_list)
        ###print("Pvalue=")
        ###print(pval_list)
        candidate_template_df <- data.frame(pval_list,distance_list,names(distance_list),i,stringsAsFactors=F)
        colnames(candidate_template_df) <- c("P_val","Distance","Template","Cell_Line")
        candidate_template_df$Q_val <- p.adjust(candidate_template_df$P_val,method="BH")
        candidate_template_df$Signif <- candidate_template_df$Q_val <= 0.1
        sub_df <- dplyr::filter(candidate_template_df, Signif==TRUE)
  ##      str(candidate_template_df)
  ##      str(sub_df)
        if(nrow(sub_df) >=1){
            ##print("Working 1")
            candidate_template_df$Prediction <- sub_df[which(sub_df$Distance == min(sub_df$Distance)),"Template"]} else if (nrow(sub_df) == 0){
                ##print("Working 2")
                candidate_template_df$Prediction <- candidate_template_df[which(candidate_template_df$Distance == min(candidate_template_df$Distance,na.rm=TRUE)),"Template"]}
##      str(candidate_template_df)

##        str(random_list)
        predicted_template_df <- rbind(predicted_template_df,candidate_template_df)                                      
                                        #str(sub_df)
        if(verbose==TRUE){
        print(paste0("Finished ", grep(paste0("^",i,"$"), colnames(alt_df)), " of ",ncol(alt_df)))} 
    }

    output <- list(templates,predicted_template_df,sample_template_list)
    names(output) <- c("Templates","Predictions","Reference_Templates")
    print("Finished")
        return(output)
}

get_template_genes_v2 <- function(expression_df, tsne_df, cluster_column="Cluster",is_ranked=TRUE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,value_returned=c("Templates","Matrix"),variable_length=TRUE,zscore_cutoff=1.8,abs_val=FALSE,sig_limit=1000){
    require(randomForest)
    require(dplyr)
    source("~/Andre_F_functions.R")

    sample_column <- grep("sample", colnames(tsne_df),ignore.case=TRUE, value=TRUE)
    tsne_df$Sample <- tsne_df[,sample_column]
    tsne_df <- tsne_df[which(tsne_df$Sample %in% colnames(expression_df)),]

    
    sub_df <- expression_df[,intersect(tsne_df$Sample, colnames(expression_df))]
    sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]

    if(method == "simple"){
        df <- zscore_df(do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample]))),"row")
##        str(df)

        df_abs <- do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample])))
        colnames(df) <- unique(tsne_df[,cluster_column])
        colnames(df_abs) <- unique(tsne_df[,cluster_column])
        if(is_ranked == TRUE){
            print("Working 1")
            if(variable_length==FALSE){
                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])} else if(variable_length==TRUE){

                    if(abs_val ==TRUE){ df <- abs(df)*-1
                                        

                                    } else{ df <- df}
                    templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]<=zscore_cutoff),x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))}

            names(templates) <- names(templates_abs)

                    
            
        } else if (is_ranked==FALSE){
            print("Working 2")
            if(variable_length==FALSE){

                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                                if(abs_val ==TRUE){ df <- abs(df)} else{ df <- df}
                                          templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))
                



                                      }
        names(templates) <- colnames(df)
        
                                  }}

    else if(method == "randomforest"){
        rownames(sub_df) <- gsub("-","_",rownames(sub_df))
        ##str(sub_df)
##        sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),] 
    transposed_df <- as.data.frame(t(sub_df))
    transposed_df[,cluster_column] <- as.factor(tsne_df[colnames(expression_df),cluster_column])
    uniq_classes <- unique(tsne_df[,cluster_column])
##    print(index)
##   print(cluster_column)


##    str(transposed_df)
    print("Starting Random Forest")
    rf <- randomForest(as.formula(paste0(cluster_column,"~.")),data=transposed_df,importance=TRUE)
        importance_df <- as.data.frame(importance(rf))[setdiff(colnames(importance(rf)),uniq_classes)]
##        str(importance_df)

    importance_df <- importance_df[order(-importance_df[,2]),]
    final_importance_matrix <- data.frame(importance(rf)[rownames(importance_df),uniq_classes],stringsAsFactors=F)

    
    colnames(final_importance_matrix) <- uniq_classes
##        str(final_importance_matrix)
##        str(colnames(final_importance_matrix))
        
    templates_abs <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, sort_order=TRUE))[1:sig_limit])

    templates_rel <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(zscore_df(final_importance_matrix),x,sort_order=TRUE))[1:sig_limit])
#    str(templates_rel)
        ##str(final_importance_matrix)
    names(templates_abs) <- colnames(final_importance_matrix)
    names(templates_rel) <- colnames(final_importance_matrix)

  ##      str(templates_abs)
  ##      str(templates_rel)
     if(variable_length==FALSE){
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                zscore_importance_matrix <- zscore_df(final_importance_matrix,"row",central_tendency="median")
                templates_rel <- lapply(colnames(zscore_importance_matrix), function(x) names(sort(zscore_importance_matrix[which(zscore_importance_matrix[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                names(templates_rel) <- colnames(zscore_importance_matrix)                
##                str(templates_rel)
##                str(templates_abs)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))      }


#    templates <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, decreasing=TRUE))[1:template_length])
    templates <- lapply(templates, function(x) gsub("_","-",x))
        names(templates) <- names(templates_abs)

    
}
#    str(templates)
        templates <- lapply(templates, function(x) x[!is.na(x)])
    if(value_returned=="Templates"){ return(templates)} else if(value_returned=="Matrix"){ return(final_importance_matrix)}
    
}




nearest_template_prediction_v3 <- function(reference_df,alt_df, tsne_df, cluster_column="Cluster",is_ranked=TRUE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,verbose=TRUE,random_sample=100,variable_length=TRUE,zscore_cutoff=1.4,templates=NULL, template_to_1=FALSE,sample_col="Sample") {
    require(lsa)
    require(dplyr)

    tsne_df <- dplyr::filter(tsne_df, get(cluster_column) !="")
    predicted_template_df <- data.frame(stringsAsFactors=F)


    template_df <- unique(tsne_df[unique(c(cluster_column))])

    template_df <- apply(template_df,2, as.character)
##    gene_set <- unique(unlist(lapply(unique(tsne_df[,TF_column]), function(x) unlist(strsplit(x,"-")))))
    gene_set <- rownames(alt_df)


    if(is.null(templates)==TRUE){
    templates <- get_template_genes_v2(reference_df, tsne_df,cluster_column,quantile_cutoff=quantile_cutoff,method=method, template_length=template_length, is_ranked=is_ranked,value_returned="Templates",variable_length=variable_length, zscore_cutoff=zscore_cutoff)} else { templates <- lapply(templates, function(x) intersect(x, rownames(alt_df))) }

    if(!is.null(reference_df)){
    if(!all(unique(unlist(templates)) %in% union(rownames(reference_df), rownames(alt_df)))){ print("There are template genes missing from the query dataframe. Dropping missing genes")
    templates <- lapply(templates, function(x) intersect(x, intersect(rownames(reference_df),rownames(alt_df))))} else{ templates <- templates}

    print("Finished templates")


###########

    samples <- lapply(sort(unique(template_df[,cluster_column])), function(x) dplyr::filter(tsne_df, get(cluster_column)==x)$Sample)
    names(samples) <- sort(unique(template_df[,cluster_column]))
    samples <- lapply(samples, function(x) intersect(x, colnames(reference_df)))



    sample_template_list <- lapply(names(samples), function(x) rowMeans(reference_df[templates[[x]],samples[[x]]],na.rm=TRUE))} else if(is.null(reference_df)){

                                                                                                                                  templates2 <- lapply(templates, function(x) unique(c(intersect(x, rownames(alt_df)), setdiff(x, x=unlist(templates)))))
                                                                                                                                  names(templates2) <- names(templates)

                                                                                                                                  sample_template_list <- list()
                                                                                                                                  for(i in names(templates)){
                                                                                                                                      vec <- templates2[[i]]
                                                                                                                                      index <- which(vec %in% templates[[i]])
                                                                                                                                      out_vec <- rep(0,length(vec)); names(out_vec) <- vec
                                                                                                                                      out_vec[index] <- 1
                                                                                                                                                                                                                                                                       sample_template_list[[i]] <- out_vec}

                                                                                                                                  templates <- templates2 }

    str(names(templates))


#######################
#######################

    for(i in colnames(alt_df)){

        print(i)

        distance_list <- unlist(lapply(sort(unique(names(templates))), function(x) 1-cosine(alt_df[templates[[x]],i],sample_template_list[[x]])))
        names(distance_list) <- sort(unique(names(templates)))


        random_list <- list()
        for(k in names(templates)){
            random_list[[k]] <- unlist(lapply(1:random_sample, function(x) 1-cosine(alt_df[sample(gene_set, length(sample_template_list[[k]])),i],sample_template_list[[k]])))
            random_list[[k]] <- random_list[[k]][!is.na(random_list[[k]])]



        }



        pval_list <- unlist(lapply(names(distance_list), function(x) 1-pnorm(distance_list[[x]],mean(random_list[[x]]),sd(random_list[[x]]),lower.tail=F)))
##        str(pval_list)
        names(pval_list) <- names(distance_list)

        candidate_template_df <- data.frame(pval_list,distance_list,names(distance_list),i,stringsAsFactors=F)
        colnames(candidate_template_df) <- c("P_val","Distance","Template","Cell_Line")
        candidate_template_df$Q_val <- p.adjust(candidate_template_df$P_val,method="BH")
        candidate_template_df$Signif <- candidate_template_df$Q_val <= 0.1
        sub_df <- dplyr::filter(candidate_template_df, Signif==TRUE)

        if(nrow(sub_df) >=1){

            candidate_template_df$Prediction <- sub_df[which(sub_df$Distance == min(sub_df$Distance)),"Template"]} else if (nrow(sub_df) == 0){

                candidate_template_df$Prediction <- candidate_template_df[which(candidate_template_df$Distance == min(candidate_template_df$Distance,na.rm=TRUE)),"Template"]}



        predicted_template_df <- rbind(predicted_template_df,candidate_template_df)                                      

        if(verbose==TRUE){
            print(paste0("Finished ", grep(paste0("^",i,"$"), colnames(alt_df)), " of ",ncol(alt_df)))}

    }

    output <- list(templates,predicted_template_df,sample_template_list)
    names(output) <- c("Templates","Predictions","Reference_Templates")
    print("Finished")
        return(output)
}


plot_key_TF_essentiality <- function(Achilles_df, NTP_results,TSNE_df){
    require(reshape2)
    require(dplyr)

    rownames(NTP_results) <- NTP_results[,1]
    TF_df <- unique(TSNE_df[c("Disease","Signif_drivers")])
    rownames(TF_df) <- TF_df[,1]

    Achilles_df <- Achilles_df[,intersect(colnames(Achilles_df),NTP_results[,1])]
    Achilles_df <- reshape2::melt(as.matrix(Achilles_df))

    Achilles_df[1:2] <- apply(Achilles_df[1:2],2,as.character)
    colnames(Achilles_df)[1:2] <- c("Gene","Cell_Line")
    Achilles_df$Cancer <- NTP_results[Achilles_df[,2],"Subtype"]
    Achilles_df$Prediction <- NTP_results[Achilles_df[,2],"Prediction"]
    Achilles_df$Signif <- NTP_results[Achilles_df[,2],"Signif"]

    Achilles_df <- Achilles_df %>% group_by(Prediction) %>% mutate(Key_TF=grepl(Gene, TF_df[Prediction,2])) %>% data.frame
    str(Achilles_df)}


sort_df <- function(df,column,sort_order=c(TRUE,FALSE)){
    df <- df[order(df[,column],decreasing=sort_order),]
    return(df)}


parse_NTP <- function(NTP_list, metadata_df,signif_cutoff=0.1,is_patient=FALSE){
    require(dplyr)

    if(is.data.frame(NTP_list)){
  
        NTP_list <- list(NTP_list)

        names(NTP_list) <- "NTP" } else{
            NTP_list <- NTP_list
    }

    for(i in names(NTP_list)){
        sub <- NTP_list[[i]]
        ##str(sub)
        sub$Signif <- sub$Q_val <= signif_cutoff
        sub_df <- dplyr::filter(sub, Signif==TRUE)
        for(k in unique(sub$Cell_Line)){
##            str(k)
            index <- which(sub$Cell_Line==k)
            int <- dplyr::filter(sub_df, Cell_Line==k)
##            str(int)
            if(nrow(int) >=1){
            ##print(int[which(int$Distance == min(int$Distance)),"Template"] == unique(sub[index,"Prediction"]))
            sub[index,"Prediction"] <- int[which(int$Distance == min(int$Distance)),"Template"]} else if (nrow(int) == 0){
#                print("Changing")
                sub[index,"Prediction"] <- sub[index,"Prediction"]




                                                                                               }}

        if(is_patient==FALSE){
        sub$Subtype <- metadata_df[sub$Cell_Line,"TCGA"]
        sub$Detailed_Cancer <- metadata_df[sub$Cell_Line,"lineage_subtype"]
        sub$Detailed_Cancer_Ext <- metadata_df[sub$Cell_Line,"disease_subtype"]
        sub$Primary <- metadata_df[sub$Cell_Line,"primary_or_metastasis"]
        sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff
        } else if(is_patient==TRUE){
            sub$Subtype <- metadata_df[sub$Cell_Line,"Subtype"]
            sub$Detailed_Cancer <- "Unknown"
            sub$Detailed_Cancer_Ext <- "Unknown"
            sub$Primary <- "Unknown"
            sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff

                                                                         } else if(is_patient =="skip"){
                                                                             sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
                                                                             sub$Q_val_cutoff <- signif_cutoff
                                                                             }

        NTP_list[[i]] <- sub
    }
    if(is_patient =="skip"){ print("Skipping accuracy assessment")
            output <- do.call("rbind",lapply(NTP_list, function(x) unique(x[c(4,6:ncol(x))])))
            final <- list(NTP_list,output)
            return(final)
            stop()} else{ 
##    str(NTP_list)
##str(NTP_list)

#print(colnames(output))

    output <- output%>% group_by(Subtype, Method, Primary) %>% mutate(Accuracy_fraction=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method) %>% mutate(Accuracy_fraction_overall=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method,Signif) %>% mutate(Accuracy_fraction_signif=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output$Primary <- gsub("^$","Unknown",output$Primary)
##output <- as.data.frame(output, stringsAsFactors=F)
##str(output)
##str(NTP_list)
final <- list(NTP_list,output)
##str(final)
    return(list(NTP_list,output))
        }}

plot_NTP_accuracy <- function(NTP_df,outfile,metadata_df,column="Disease",plot_width=24,plot_height=8,ncol=5){
    require(dplyr)

    print("Finished with the easy stuff")
    #str(NTP_df)
    uniq_methods <- unique(NTP_df$Method)
    print(uniq_methods)
  
    sub <- lapply(unique(NTP_df$Method), function(x) table(dplyr::filter(NTP_df, Method==x)[c("Prediction","Subtype")]))
    #str(sub)
    sub_summary <- sub
                                        #    str(NTP_df)
                           
    sub_signif <- lapply(unique(NTP_df$Method), function(x) table(dplyr::filter(NTP_df, Method==x)[c("Prediction","Subtype","Signif")]))
##    str(sub_signif)
  
    for(i in 1:length(sub)){ int <- as.data.frame(rowSums(sub[[i]]))
                             colnames(int) <- "Freq"
                             int$Template <- rownames(int)
                             int$Method <- unique(NTP_df$Method)[i]
                             sub_summary[[i]] <- int

                             int2 <- as.data.frame(sub[[i]])
                             ##str(int2)
                             colnames(int2) <- c("Template","Cancer","Freq")
                             
                             int2$Method <- unique(NTP_df$Method)[i]
                             
                             sub[[i]] <- int2

                             int3 <- as.data.frame(sub_signif[[i]],stringsAsFactors=F)
                             
                             colnames(int3) <- c("Template","Cancer","Signif","Freq")
                             int3$Method <- unique(NTP_df$Method)[i]
                             ##str(int3)
                             sub_signif[[i]] <- int3
                         }

    sub_summary <- do.call("rbind",sub_summary)
    sub <- do.call("rbind",sub)
    sub_signif <- do.call("rbind",sub_signif)
    ##str(sub_summary)
    ##str(sub)
    ##print(summary(sub$Freq))
    ##str(sub)

    p <- ggplot(NTP_df, aes(Subtype, Accuracy_fraction_overall))+geom_col(position="dodge")+facet_grid(Method~.)+theme_bw()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=12))+geom_hline(yintercept=0.6, color="gray50", linetype=2,size=2)
    
    p1 <- ggplot(NTP_df, aes(Subtype, Accuracy_fraction,fill=Primary))+geom_col(position="dodge")+facet_grid(Method~.)+theme_bw()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=12))+geom_hline(yintercept=0.6, color="gray50", linetype=2,size=2)
    
    
    p2 <- ggplot(NTP_df, aes(Method, Accuracy_fraction_overall,fill=as.factor(Method)))+geom_col(position="dodge")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=15,angle=30,hjust=1),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=12))+geom_hline(yintercept=0.6, color="gray50", linetype=2,size=2)+facet_wrap(Subtype~.,ncol=ncol)+labs(title="Accuracy of all predictions")
    p2_5 <- ggplot(dplyr::filter(NTP_df,Signif==TRUE), aes(Method, Accuracy_fraction_signif,fill=as.factor(Method)))+geom_col(position="dodge")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=15,angle=30,hjust=1),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=12))+geom_hline(yintercept=0.6, color="gray50", linetype=2,size=2)+facet_wrap(Subtype~.,ncol=ncol)+labs(title="Accuracy of significant predictions")
    

    p3 <- ggplot(NTP_df,aes(Subtype, Prediction, label=Cell_Line, color=Subtype,alpha=0.7))+geom_jitter(size=3,width=0.1,height=0.15)+geom_text(fontface="bold",size=3,check_overlap=TRUE,show.legend=F)+facet_wrap(Method~.,ncol=2)+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10),axis.text.y=element_text(face='bold',size=10),strip.text = element_text(colour = "black", face = "bold",size=12))+scale_TF_score("Corces_Hex_Codes_Cancer.txt")
    

    p4 <- ggplot(NTP_df,aes(Subtype, Prediction, label=Cell_Line, color=Subtype,shape=Primary,alpha=0.7))+geom_jitter(size=3,width=0.1,height=0.15)+geom_text(fontface="bold",size=3,check_overlap=TRUE,show.legend=F)+facet_wrap(Method~.,ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10),axis.text.y=element_text(face='bold',size=10),strip.text = element_text(colour = "black", face = "bold",size=8))+scale_TF_score("Corces_Hex_Codes_Cancer.txt")
    
    p5 <- ggplot(NTP_df, aes(Prediction,fill=Subtype,label=Subtype,alpha=0.8))+geom_histogram(stat="count")+scale_TF_score("Corces_Hex_Codes_Cancer.txt")+facet_grid("Method")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=16))+facet_wrap(Method~.,ncol=1)+labs(title="Template assignments: All")

    p5_5 <- ggplot(dplyr::filter(NTP_df,Signif==TRUE), aes(Prediction,fill=Subtype,label=Subtype,alpha=0.8))+geom_histogram(stat="count")+scale_TF_score("Corces_Hex_Codes_Cancer.txt")+facet_grid("Method")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=16))+facet_wrap(Method~.,ncol=1)+labs(title="Template assignments: Significant")

    p5_75 <- ggplot(dplyr::filter(NTP_df,Signif==FALSE), aes(Prediction,fill=Subtype,label=Subtype,alpha=0.8))+geom_histogram(stat="count")+scale_TF_score("Corces_Hex_Codes_Cancer.txt")+facet_grid("Method")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=16))+facet_wrap(Method~.,ncol=1)+labs(title="Template assignments: **NOT** Statistically Significant")

print("Working Heatmaps")

#####################################################################################################################
    
    centroid_df <- as.data.frame(do.call("rbind",lapply(unique(metadata_df[,column]), function(x) colMeans(dplyr::filter(metadata_df, get(column)==x)[,1:2]))),stringsAsFactors=F)
    centroid_df[,column] <- unique(metadata_df[,column])

    centroid_df <- centroid_df[hclust(dist(centroid_df[,1:2]))$order,]
    factor_levels <- centroid_df[,column]
    sub$Template <- factor(sub$Template, levels=factor_levels)
    sub_signif$Template <- factor(sub_signif$Template, levels=factor_levels)

    centroid_df2 <- tryCatch(as.data.frame(do.call("rbind",lapply(unique(metadata_df[,"Subtype"]), function(x) colMeans(dplyr::filter(metadata_df, Subtype==x)[,1:2]))),stringsAsFactors=F))

    if(is.data.frame(centroid_df2) ==TRUE && nrow(centroid_df2) >=2){
        centroid_df2$Subtype <- unique(metadata_df$Subtype)
        centroid_df2$ordering <- unlist(lapply(1:nrow(centroid_df2),function(x) min(grep(centroid_df2[x,3], centroid_df[,3]))))
        centroid_df2 <- sort_df(centroid_df2,"ordering",sort_order=FALSE)
        factor_levels2 <- centroid_df2[,3]

        ###Quick hack of factor levels for better plots #####
        factor_levels2 <- c("LIHC", "TGCT","LUAD","LUSC","STAD","ESCA","CHOL","THCA","COAD","BRCA","SKCM","BLCA","UCEC","HNSC","CESC","PRAD","KIRC","ACC","MESO","PCPG","LGG")
       
        sub$Cancer <- factor(as.character(sub$Cancer), levels=factor_levels2)
        sub_signif$Cancer <- factor(as.character(sub_signif$Cancer), levels=factor_levels2)
        print("Finished arranging heatmaps")
    } else{print("Ordering columns did not work")}


    sub <- sub %>% group_by(Cancer,Method) %>% mutate(Freq_Norm=round(Freq/sum(Freq),digits=2)) %>% data.frame
    sub_signif <- sub_signif %>% group_by(Cancer,Method) %>% mutate(Freq_Norm=round(Freq/sum(Freq),digits=2)) %>% data.frame

#####################################################################################################################    
    p7 <- ggplot(NTP_df, aes(Num_Predictions,fill=as.factor(Method),alpha=0.4))+geom_histogram(binwidth=1,position="identity")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15),strip.text = element_text(colour = "black", face = "bold",size=16))+facet_grid(Method~.)

    p8 <- ggplot(sub, aes(y=Template,Cancer,fill=Freq,label=Freq))+geom_tile()+geom_text(fontface="bold",size=6,color="gray50")+scale_fill_viridis_c(option="magma")+facet_wrap(Method~., ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="All Predictions: Raw Numbers")
    
    p9 <- ggplot(sub, aes(y=Template,Cancer,fill=Freq_Norm,label=Freq_Norm))+geom_tile()+geom_text(fontface="bold",size=4,color="gray50")+scale_fill_viridis_c(option="magma")+facet_wrap(Method~., ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="All Predictions: Normalized Numbers")

    print("Almost done")
     
   p10 <- ggplot(dplyr::filter(sub_signif,Signif==TRUE), aes(y=Template,Cancer,fill=Freq,label=Freq))+geom_tile()+geom_text(fontface="bold",size=6,color="gray50")+scale_fill_viridis_c(option="magma")+facet_wrap(Method~., ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Raw Numbers")
    
    sub_signif2 <- dplyr::filter(sub_signif,Signif==TRUE)
    sub_signif2 <- sub_signif2 %>% group_by(Cancer,Method) %>% mutate(Freq_Norm=round(Freq/sum(Freq),digits=2)) %>% data.frame
##    str(sub_signif2)
    print("One.. Last.. Plot")
##    print(summary(sub_signif2[c("Freq_Norm","Method")]))
    p11 <- ggplot(sub_signif2, aes(y=Template,Cancer,fill=Freq_Norm,label=Freq_Norm))+geom_tile()+geom_text(fontface="bold",size=4,color="gray50")+scale_fill_viridis_c(option="magma")+facet_wrap(Method~., ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")

    p12 <- ggplot(sub_signif2, aes(y=Template,Cancer,fill=Freq_Norm,label=Freq))+geom_tile()+geom_text(fontface="bold",size=4,color="gray50")+scale_fill_viridis_c(option="magma")+facet_wrap(Method~., ncol=2,scales="free")+theme_bw()+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")

    pdf(outfile,width=plot_width,height=plot_height)

    print(p)
    print(p1)
    print(p2)
    print(p2_5)
    print(p3)
    print(p4)
    print(p5)
    print(p5_5)
    try(print(p5_75),silent=TRUE)
    print(p7)
    print(p8)
    print(p9)
    print(p10)
    print(p11)
    print(p12)
    dev.off()
    print("All Done")

}

make_coessential_corplot <- function(Achilles_df, gene_set, cor_method=c("pearson","spearman")){
    require(dplyr)
    require(ggplot2)
    
    combinations <- lapply(1:length(gene_set)^2,function(x) unlist(expand.grid(gene_set, gene_set, stringsAsFactors=F)[x,]))

    df<- do.call("rbind",lapply(combinations, function(x) data.frame(unlist(Achilles_df[x[1],]), unlist(Achilles_df[x[2],]), x[1],x[2],stringsAsFactors=FALSE,colnames(Achilles_df))))
    colnames(df) <- c("Essentiality1","Essentiality2","Gene1","Gene2","Sample")
    df <- df %>% group_by(Gene1,Gene2) %>% mutate(Corr=round(suppressWarnings(cor.test(Essentiality1,Essentiality2,method=cor_method))$estimate,digits=3), P_val=round(suppressWarnings(cor.test(Essentiality1,Essentiality2,method=cor_method))$p.value,digits=3)) %>% data.frame
    df$Corr_Method <- cor_method
    str(df)
    df$Corr_Method <- gsub("spearman","Spearman",df$Corr_Method)
    df$Corr_Method <- gsub("pearson", "Pearson",df$Corr_Method)
    df$Facet <- paste0(df$Gene1,"~",df$Gene2)
    df$Signif <- df$P_val <=0.05
    df$Flag <- df$Gene1 == df$Gene2
##    df <- dplyr::filter(df, Flag==FALSE)


    anno_df <- unique(df[setdiff(colnames(df), c("Essentiality1","Essentiality2","Sample"))])
    anno_df <- anno_df %>% mutate(Anno=paste0("Corr:",Corr,"\nP-value:",P_val))
##    str(anno_df)

    p <- ggplot(df,aes(Essentiality1,Essentiality2,label=Sample))+geom_text(fontface="bold", size=0.3)+geom_point(aes(color=Signif),shape=19,size=2,show.legend=F)+geom_smooth(method="lm",color="gray50",show.legend=F)+xlab(NULL)+ylab(NULL)+theme_classic()+theme(axis.text.x=element_text(face='bold',size=10),axis.text.y=element_text(face='bold',size=10),strip.background=element_rect(color="black",fill="gray"),strip.text = element_text(colour = "black", face = "bold",size=10))+facet_grid(Gene1~Gene2, scales="free")+geom_text(data=anno_df, aes(x=-Inf,y=Inf,vjust=0,hjust=1,label=Anno,alpha=1),fontface="bold",size=3,show.legend=F)+labs(title=paste0(unique(anno_df$Corr_Method)," correlation across samples"))+scale_color_manual(values=c("black","red","white"))+scale_y_reverse()+scale_x_reverse()
    return(p)


}
    
get_coessential_genes_v2 <- function(Achilles_df, metadata_df=NULL,subset=NULL, genes_of_interest=rownames(Achilles_df),precomputed_mat=NULL,corr_cutoff=0.9,cor_method=c("pearson","spearman")){
    require(Hmisc)
    require(dplyr)
    require(reshape2)
    print(dim(Achilles_df))
    genes_of_interest <- intersect(genes_of_interest, rownames(Achilles_df))
   
    if(is.null(metadata_df) && is.null(subset)){
        print("Working with all data")
        Achilles_df <- Achilles_df
        flag <- "All"

    } else {
            Achilles_df <- Achilles_df[,intersect(colnames(Achilles_df),metadata_df[,1])]
            metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(Achilles_df)),]
            index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
            print(subset)
            print(index)
            print(metadata_df[index,1])

            ##stop()

            if(length(index) <5 & length(index)>3){ print(paste0("Not enough observations for very reliable analysis, proceeding with subset + ",5-length(index)," random samples"))
                                                    index <- c(index, sample(setdiff(1:nrow(metadata_df), index),5-length(index)))
                                                    flag <- "Subset"
                                                    Achilles_df <- Achilles_df[metadata_df[index,1]]

                                                } else if (length(index) <= 3){ print("0 to 3 cell lines match your search criteria, Running in global correlation mode")
                                                                                ##stop()
                                                                                Achilles_df <- Achilles_df
                                                                                flag <- "Subset"
                                                                } else{
                                                                    Achilles_df <- Achilles_df[metadata_df[index,1]]
    flag <- "Subset"
                                                                }}
    if(is.null(precomputed_mat)==TRUE){
        print("Calculating correlation matrix")
        cormat <- rcorr(t(Achilles_df),type=cor_method)} else if(is.null(precomputed_mat)==FALSE){ print("Working with precomputed matrix")
                                                cormat <- precomputed_mat}
    print("Finished Correlation matrix")

    gc()
    
    cor_vals <- reshape2::melt(cormat[[1]][genes_of_interest,])
    cor_vals[,1:2] <- apply(cor_vals[,1:2],2,as.character)
    print("Finished melting correlation matrix")

    cor_pvals <- reshape2::melt(cormat[[3]][genes_of_interest,])
    cor_pvals[,1:2] <- apply(cor_pvals[,1:2],2,as.character)
    print("Finished melting p-value matrix")
#        str(cor_vals)
#        str(cor_pvals)
        cor_vals$pvals <- cor_pvals[,3]

        sub_vals <- dplyr::filter(cor_vals,pvals <=0.1)
        out_list <- list()
    print("Finished melting")
    str(genes_of_interest)
        for(i in genes_of_interest){
            sub <- dplyr::filter(sub_vals, Var1==i)
            ##sub <- sub_vals
            sub$qvals <- p.adjust(sub$pvals,method="BH")
            self <- data.frame(cbind(sub[1,1],sub[1,1],1,0,0))
            colnames(self) <- colnames(sub)
#            str(self)
#           str(sub)
            sub <- rbind(sub,self,stringsAsFactors=F)
            sub[3:5] <- apply(sub[3:5],2,as.numeric)
            out_list[[i]] <- dplyr::filter(sub, qvals <=0.2)
            print(paste0("Finished coessentiality analysis for ",i))
        }

        out_list <- unique(do.call("rbind",out_list))
    str(out_list)
    colnames(out_list) <- c("Gene1","Gene2","Correlation","P_value","Q_value")
    out_list[,1] <- as.character(out_list[,1]); out_list[,2] <- as.character(out_list[,2])
    out_list$Quantile <- ecdf(abs(out_list$Correlation))(abs(out_list$Correlation))
    out <- dplyr::filter(out_list, Quantile>corr_cutoff)
    print("All done with coessentiality")

    return(list(cormat,out_list,out,length(index),flag))}
                          
get_final_candidate_scores_v2 <- function(druggable_genes_output, TF_score_output, essentiality_df, CGC_list, TF_list,flag=c("All","Subset"),subset=NULL,cell_line_mapping_df=NULL){
    require(dplyr)

    if(is.null(cell_line_mapping_df) ==FALSE){ index <- intersect(dplyr::filter(cell_line_mapping_df, Prediction==subset,Signif==TRUE)[,1],colnames(essentiality_df))
                                               if(length(index)<1){
                                                   index <- colnames(essentiality_df)} else {index <- index}
                                               str(index)

                                           } else{ print("Not calculating within cluster essentiality")
                                   index <- colnames(essentiality_df)                                                                                                                                                          }
#######Getting essentiality scores for Genes, obtaining p-value of median essentiality relative to others#########
    druggable_genes_output$Essentiality <- unlist(lapply(druggable_genes_output$Gene2, function(x) median(unlist(essentiality_df[x,]),na.rm=TRUE)))
    if(length(index)!=ncol(essentiality_df)){
        druggable_genes_output$Essentiality_in_cluster <- unlist(lapply(druggable_genes_output$Gene2, function(x) median(unlist(essentiality_df[x,index]),na.rm=TRUE)))} else {druggable_genes_output$Essentiality_in_cluster <- 0}
    
    druggable_genes_output$Essentiality_pval <- 2*(1-pnorm(abs(druggable_genes_output$Essentiality),mean(druggable_genes_output$Essentiality,na.rm=TRUE),sd(druggable_genes_output$Essentiality,na.rm=TRUE)))
    druggable_genes_output$Always_essential <- ifelse(druggable_genes_output$Essentiality <0 & druggable_genes_output$Essentiality_pval <0.05,TRUE,FALSE)
print("Finished Essentiality")
########Get TF score clusters for each gene #####################
    druggable_genes_output$Cluster_Gene1 <- unlist(lapply(druggable_genes_output$Gene1,function(x) paste0(unique(TF_score_output[grep(paste0("^",x,"-|-",x,"-|-",x,"$"), TF_score_output$Signif_drivers),"Disease"]),collapse="|")))
    druggable_genes_output$Cluster_Gene2 <- unlist(lapply(druggable_genes_output$Gene2,function(x) paste0(unique(TF_score_output[grep(paste0("^",x,"-|-",x,"-|-",x,"$"), TF_score_output$Signif_drivers),"Disease"]),collapse="|")))

    druggable_genes_output$Cluster_TF1 <- unlist(lapply(druggable_genes_output$Cluster_Gene1, function(x) paste0(unique(dplyr::filter(TF_score_output, Disease %in% unlist(strsplit(x,"\\|")))$Signif_drivers),collapse="|")))
    druggable_genes_output$Cluster_TF2 <- unlist(lapply(druggable_genes_output$Cluster_Gene2,function(x) paste0(unique(dplyr::filter(TF_score_output, Disease %in% unlist(strsplit(x,"\\|")))$Signif_drivers),collapse="|")))

    print("Finished getting Clusters")
    
####### Get number of drugs gene recurrence and number of clusters per gene ########
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Num_Clusters=length(unique(unlist(strsplit(Cluster_Gene1, "\\|"))))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Drugs) %>% mutate(Num_Drugs=length(unique(unlist(strsplit(Drugs,"~"))))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene2) %>% mutate(Gene_recurrence=length(unique(Gene1))) %>% data.frame
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Targetable=(Num_Drugs>=1)) %>% data.frame

############### Check if Gene is in CGC list or a TF ########
    druggable_genes_output$Is_TF <- druggable_genes_output$Gene2 %in% TF_list
    druggable_genes_output$CGC <- druggable_genes_output$Gene2 %in% CGC_list
   
########### Get Mean rank for TF in clusters where it's important #####
    print(paste0("Getting TF score ranks for ", subset,"... This could take a while"))
##    saveRDS(druggable_genes_output,"test_hlist.rds")
    druggable_genes_output$TF_score_rank <- make_hierarchical_list(druggable_genes_output,flag,subset,search_index=1)
    str(druggable_genes_output)
    druggable_genes_output$Gene_TF_score_rank <- make_hierarchical_list(druggable_genes_output,flag,subset,search_index=2)
    print("Finalizing scores")
    
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Simple_score=Targetable*(1*(abs(Correlation))-abs(Essentiality)+(1/Num_Clusters))) %>% mutate(Simple_score_with_TF_rank=Targetable*Simple_score*((1/TF_score_rank))) %>% mutate(Simple_score_TF_rank_Tox=Targetable*Simple_score*(1/TF_score_rank)*(1-Toxic_fraction)) %>% mutate(Simple_score_TF_rank_is_TF_is_Tox=Targetable*Simple_score*((Is_TF+1-Toxic_fraction)*((1/TF_score_rank)))) %>% data.frame

    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(New_Simple_score=Targetable*(1*(abs(Correlation))-abs(Essentiality)-Essentiality_in_cluster)*(1+1/Num_Clusters)) %>% mutate(New_Simple_score_with_TF_rank=Targetable*New_Simple_score*(1+(1/TF_score_rank))) %>% mutate(New_Simple_score_TF_rank_Tox=Targetable*New_Simple_score*(1+(1/TF_score_rank))*(1-Toxic_fraction)) %>% mutate(New_Simple_score_TF_rank_is_TF_is_Tox=Targetable*New_Simple_score*(1-Toxic_fraction)*(Is_TF+(1/TF_score_rank))) %>% mutate(Optimized_Score=Targetable*(1*(abs(Correlation))-abs(Essentiality)-Essentiality_in_cluster-((Essentiality_in_cluster-Essentiality)))*(1+1/Num_Clusters)*(1+(1/TF_score_rank))*(1-Toxic_fraction)) %>% mutate(Optimized_Score_No_Tox=Targetable*(1*(abs(Correlation))-abs(Essentiality)-Essentiality_in_cluster-((Essentiality_in_cluster-Essentiality)))*(1+1/Num_Clusters)*(1+(1/TF_score_rank))) %>% mutate(Optimized_Score_Tweak_Tox=Targetable*(1*(abs(Correlation))-abs(Essentiality)-Essentiality_in_cluster-((Essentiality_in_cluster-Essentiality)))*(1+1/Num_Clusters)*(1+(1/TF_score_rank))*(1+(1-Toxic_fraction))) %>% data.frame
    
    druggable_genes_output[is.na(druggable_genes_output)] <- 0

    druggable_genes_output <- sort_df(druggable_genes_output,"Optimized_Score",sort_order=TRUE)
    druggable_genes_output <- dplyr::filter(druggable_genes_output, TF_score_rank <1000)
    druggable_genes_output$Cancer_Type <- subset
    druggable_genes_output$Num_Cell_Lines <- length(index)
    druggable_genes_output$Mode <- flag
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene2) %>% mutate(Candidates=paste0(Cancer_Type,"~",Gene2),Final_Score=max(Optimized_Score)) %>% ungroup %>% data.frame

    rank_df <- unique(druggable_genes_output[c("Gene2","Final_Score")])
    rank_df$Rank <- (1+nrow(rank_df))-rank(unlist(rank_df[,2]))
    rownames(rank_df) <- rank_df[,1]
    druggable_genes_output$Final_Rank <- rank_df[druggable_genes_output$Gene2,"Rank"]
    druggable_genes_output$Top3 <- paste0(unique(dplyr::filter(druggable_genes_output,Final_Rank<=3)$Gene2),collapse="-")
    
    print("All Done")

return(druggable_genes_output)

}


make_hierarchical_list <- function(df,flag=c("All","Subset"),subset,column1="Cluster_TF1",column2="Cluster_Gene1",search_index=c(1,2)){
    print(flag)
    print(subset)
    if(flag =="All" && is.null(subset) ==TRUE){
        print("Working 1")
        output <- unlist(lapply(1:nrow(df), function(x) mean(unlist(lapply(strsplit(df[x,"Cluster_TF1"],"\\|"), function(y) lapply(strsplit(y,"-"), function(z) grep(df[x,1],z)))),na.rm=TRUE))) }    else if(flag=="Subset" && is.null(subset)==FALSE){ 
    ##else if(flag=="All" && is.null(subset)==FALSE){     ##    print("Working 2")}
    print("Working 3")
    hlist <- list()

    for(i in 1:nrow(df)){
        ##print(i)

        pair <- paste0(df[i,1],"~",df[i,2])
        gene1 <- df[i,1]
        ##print(pair)
        

        

        int_l1 <- unlist(strsplit(df[i,"Cluster_TF1"],"\\|"))
        ##str(int_l1)
        int_l1 <- lapply(int_l1, function(x) unlist(strsplit(x,"-")))
        
        
        names_int <- unlist(strsplit(df[i,"Cluster_Gene1"],"\\|"))
        names(int_l1) <- names_int
        ##print(int_l1)
        ##print(names_int)
        
        hlist[[pair]] <- int_l1
        ##print(hlist)
    }
##    str(hlist)

    ##print(hlist)

    score_vec <- lapply(1:nrow(df), function(x) grep(df[x,search_index], hlist[[paste0(df[x,1],"~",df[x,2])]][[subset]]))
##    print(table(unlist(score_vec)))
    for(k in 1:length(score_vec)){
        val <- ifelse(length(score_vec[[k]])==0, 1000,score_vec[[k]])
        score_vec[[k]] <- val
        output <- unlist(score_vec)}}

    return(output)}

rescore_druggable_genes <- function(druggable_genes_output,subset,flag="Subset",adjust_TF_rank=TRUE){
    require(dplyr)

    if(adjust_TF_rank==TRUE){
        druggable_genes_output$TF_score_rank <- make_hierarchical_list(druggable_genes_output,flag,subset)
        druggable_genes_output$Gene_TF_score_rank <- make_hierarchical_list(druggable_genes_output,flag,subset,search_index=2)
        


    } else if(adjust_TF_rank==FALSE){ print("Not calculating new TF ranks")}
    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Num_Clusters=length(unique(unlist(strsplit(Cluster_Gene1, "\\|"))))) %>% data.frame

    druggable_genes_output <- druggable_genes_output[setdiff(1:ncol(druggable_genes_output), grep("Simple_score", colnames(druggable_genes_output)))]
##    str(druggable_genes_output)
    
        druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Simple_score=Targetable*(1*(abs(Correlation))-abs(Essentiality))*(1+(1/Num_Clusters))) %>% mutate(Simple_score_with_TF_rank=Targetable*Simple_score*((1/TF_score_rank))) %>% mutate(Simple_score_TF_rank_Tox=Targetable*Simple_score*(1/TF_score_rank)*(1-Toxic_fraction)) %>% mutate(Simple_score_TF_rank_is_TF_is_Tox=Targetable*Simple_score*(1-Toxic_fraction)*(Is_TF+(1/TF_score_rank))) %>% data.frame

    druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(New_Simple_score=Targetable*(1*(abs(Correlation))-abs(Essentiality)-Essentiality_in_cluster)*(1+(1/Num_Clusters))) %>% mutate(New_Simple_score_with_TF_rank=Targetable*New_Simple_score*((1/TF_score_rank))) %>% mutate(New_Simple_score_TF_rank_Tox=Targetable*New_Simple_score*(1/TF_score_rank)*(1-Toxic_fraction)) %>% mutate(New_Simple_score_TF_rank_is_TF_is_Tox=Targetable*New_Simple_score*(1-Toxic_fraction)*(Is_TF+(1/TF_score_rank))) %>% data.frame

   
  ## druggable_genes_output <- druggable_genes_output %>% group_by(Gene1,Gene2) %>% mutate(Simple_score=Targetable*(1*(abs(Correlation))+Essentiality+(1/Num_Clusters))) %>% mutate(Simple_score_with_TF_rank=Targetable*Simple_score*(1+(1/TF_score_rank))) %>% mutate(Simple_score_TF_rank_is_TF=Targetable*Simple_score*(Is_TF*(1+(1/TF_score_rank)))) %>% mutate(Simple_score_TF_rank_is_TF_is_Tox=Targetable*Simple_score*((Is_TF+1-Toxic_fraction)*(1+(1/TF_score_rank)))) %>% mutate(Simple_score_TF_rank_is_TF_is_CGC=Targetable*Simple_score*((Is_TF+CGC)*(1+(1/TF_score_rank)))) %>% mutate(Simple_score_TF_rank_is_TF_is_CGC_Tox=Targetable*Simple_score*((Is_TF+1-Toxic_fraction+CGC)*(1+(1/TF_score_rank)))) %>% data.frame

    druggable_genes_output[is.na(druggable_genes_output)] <- 0

    druggable_genes_output <- sort_df(druggable_genes_output,"Simple_score_with_TF_rank",sort_order=TRUE)

        print("All Done")

        return(druggable_genes_output)}


prioritize_drug_candidates <- function(druggable_genes_output,Achilles_df ,DGIDB,num_drugs=3){
    require(dplyr)

    DGIDB$interaction_types <- ifelse(nchar(DGIDB$interaction_types)<1,"Unknown",DGIDB$interaction_types)
    DGIDB <- dplyr::filter(DGIDB, interaction_types !="Unknown")


    evidence_dict <- list("TdgClinicalTrial"=1,"OncoKB"=1,"FDA"=1,"CKB"=2,"TTD"=3,"ChemblInteractions"=3,"CGI"=3,"ClearityFoundationClinicalTrial"=3,"CancerCommons"=3,"MyCancerGenome"=3,"TEND"=3,"TALC"=3,"CIViC"=3,"MyCancerGenomeClinicalTrial"=3,"GuideToPharmacologyInteractions"=3,"DoCM"=3,"ClearityFoundationBiomarkers"=4,"TEND"=4,"NCI"=4,"empty"=5)

    interaction_dict <- list("inhibitor"=5,"agonist"=-5,"antagonist"=5,"antagonist,inhibitor"=5,"antibody"=4,"agonist,antagonist"=4, "inhibitory allosteric modulator"=4,"suppressor"=4,"gating inhibitor"=3,"activator"=-3,"allosteric modulator,antagonist"=3,"agonist,allosteric modulator"=-3,"blocker"=3,"channel blocker,gating inhibitor"=3,"negative modulator"=3,"channel blocker"=3,"partial agonist"=-2,"modulator"=2,"allosteric modulator"=2,"inverse agonist"=2,"activator,channel blocker"=-2,"antisense"=2,"antisense oligonucleotide"=2,"cofactor"=2,"binder"=2 ,"inducer"=-2,"stimulator"=-2,"positive allosteric modulator"=-2,"activator,antagonist"=2,"Unknown"=1,"vaccine"=2)

    print(setdiff(names(evidence_dict), unique(DGIDB$interaction_claim_source)))
    
    DGIDB <- DGIDB %>% group_by(gene_name) %>% mutate(gene_essentiality=median(unlist(Achilles_df[gene_name,]),na.rm=TRUE)) %>% ungroup %>% group_by(drug_name) %>% mutate(Num_Genes=length(unique(gene_name)),Drug_Toxicity=median(gene_essentiality,na.rm=TRUE)) %>% ungroup %>% group_by(gene_name,drug_name) %>% mutate(Interaction_Tier=mean(unlist(interaction_dict[interaction_types])),interaction_claim_source=ifelse(nchar(interaction_claim_source)==0,"empty",interaction_claim_source), Evidence_Tier=mean(unlist(evidence_dict[interaction_claim_source]))) %>% data.frame

    print("Grouping 1")

    DGIDB <- DGIDB %>% group_by(drug_name) %>% mutate(Final_drug_score=Interaction_Tier+1-Evidence_Tier+(1-(Num_Genes/median(unique(DGIDB[c("drug_name","Num_Genes")])$Num_Genes)))+Drug_Toxicity, Final_drug_score_redux= Interaction_Tier*Evidence_Tier*(1-(Num_Genes/median(unique(DGIDB[c("drug_name","Num_Genes")])$Num_Genes)))) %>% data.frame
    print("Grouping 2")
    
    str(DGIDB)
    ##    print(table(DGIDB$gene_essentiality==DGIDB$Drug_Toxicity))
    druggable_genes_list <- lapply(unique(druggable_genes_output$Gene2), function(x) DGIDB[which(DGIDB[,1] ==x),])
    names(druggable_genes_list) <- unique(druggable_genes_output$Gene2)
    str(names(druggable_genes_list))
    druggable_genes_list <- lapply(druggable_genes_list, function(x) sort_df(x,"Final_drug_score",TRUE))

    print("Getting ranks of top drugs")

    druggable_genes_output$Priority_Drugs <- unlist(lapply(druggable_genes_output$Gene2, function(x) paste0(druggable_genes_list[[x]][,"drug_name"][1:min(c(num_drugs,length(druggable_genes_list[[x]][,"drug_name"])))],collapse="~")))
    print("Getting scores of top drugs")
    
    druggable_genes_output$Drug_Scores <- unlist(lapply(druggable_genes_output$Gene2, function(x) paste0(round(druggable_genes_list[[x]][,"Final_drug_score"][1:min(c(num_drugs,length(druggable_genes_list[[x]][,"Final_drug_score"])))],digits=3),collapse="~")))
    ##druggable_genes_output <- druggable_genes_output %>% group_by(Gene2) %>% mutate(Drugs=ifelse(Gene2 %in% names(druggable_genes_list) ,"working","None")) %>% data.frame
##paste0(druggable_genes_list[[Gene2]]$drug_name,"~")
##&& nrow(druggable_genes_list[[Gene2]])>=1
    ##print(table(druggable_genes_output$Drugs))

    return(list(druggable_genes_output,DGIDB,druggable_genes_list))
    

}


plot_drug_candidates <- function(df, highlight_list,subset=NULL,highlight_column=c("Pair","Candidates"),scale="~/ANF_Final_Cluster_key.txt"){
    require(dplyr)
    require(ggrepel)

    if( is.null(subset)== FALSE){ df <- dplyr::filter(df, grepl(subset, Disease)) } else {
        print("Not subsetting, plotting all data")
        df <- df}
    feature_df <- df[,c("Pair","Correlation","Targetable","Toxic_fraction","Essentiality","Essentiality_in_cluster","Num_Clusters","Num_Drugs","Is_TF","CGC", "TF_score_rank","Gene_TF_score_rank","Cancer_Type")]

    score_df <- df[,c("Pair","Cancer_Type",grep("score",colnames(df),value=TRUE))]
    
    feature_df <- reshape2::melt(feature_df,id.vars=c("Pair","Cancer_Type"))
    score_df <- reshape2::melt(score_df,id.vars=c("Pair","Cancer_Type"))
    colnames(feature_df)[3] <- "Feature"
    colnames(score_df)[3:4] <- c("Schema","Score")
    


    interesting_subset_feature <- dplyr::filter(feature_df, get(highlight_column) %in% highlight_list)
    interesting_subset_score <- dplyr::filter(score_df, get(highlight_column) %in% highlight_list)

##    str(feature_df)
##   str(score_df)
##    str(interesting_subset_feature)
##    str(interesting_subset_score)

   
    p <- ggplot(score_df, aes(Score, fill=Cancer_Type))+geom_histogram(position="identity")+geom_vline(data=interesting_subset_score, aes(xintercept=Score, color=Cancer_Type),size=0.4,linetype="dashed")+facet_grid(Cancer_Type~Schema, scales="free")+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=12),strip.text.y = element_text(colour = "black", face = "bold",size=8))+scale_TF_score(scale,extension=".txt")+geom_text(data=interesting_subset_score, aes(x=Score,color=Cancer_Type ,y=+Inf,label=Pair,hjust=1, angle=90, fontface="bold"),check_overlap=FALSE)

    p2 <- ggplot(feature_df, aes(value, fill=Cancer_Type))+geom_histogram(position="identity")+geom_vline(data=interesting_subset_feature, aes(xintercept=value, color=Cancer_Type),size=0.4,linetype="dashed")+facet_grid(Cancer_Type~Feature,scales="free")+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=12),axis.text.y=element_text(face='bold',size=15),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=15),strip.text.y = element_text(colour = "black", face = "bold",size=8))+scale_TF_score(scale,extension=".txt")+geom_text(data=interesting_subset_feature, aes(x=value,color=Cancer_Type ,y=+Inf,label=Pair,hjust=1, angle=90, fontface="bold"),check_overlap=FALSE)

    print(p)
    print(p2)

#return(interesting_subset_feature)
}

plot_drug_candidates_v2 <- function(df, highlight_list,subset=NULL,highlight_column=c("Pair","Candidates"),scale="~/ANF_Final_Cluster_key.txt",score_column="Optimized_Score",fontsize=3){
    require(dplyr)
    require(reshape2)

    if( is.null(subset)== FALSE){ df <- dplyr::filter(df, Cancer_Type %in% subset) } else {
        print("Not subsetting, plotting all data")
        df <- df}
    print(nrow(df))
    df <- df[-grep("^Simple_score|_is_TF_",colnames(df))]
    df <- df[which(df$Final_Score == df[,score_column]),]
    print(nrow(df))
    df$Candidates <- paste0(df$Cancer_Type,"~",df$Gene2)
    df$Pair <- paste0(df[,1],"~",df[,2])
    feature_df <- df[,c("Pair","Candidates","Correlation","Toxic_fraction","Essentiality","Essentiality_in_cluster","Num_Clusters","TF_score_rank","Cancer_Type")]

    score_df <- df[,c("Pair","Candidates","Gene2","Cancer_Type",grep("score",colnames(df),value=TRUE,ignore.case=TRUE))]
    score_df <- score_df[-grep("TF_score",colnames(score_df))]
    final_score_df <- unique(df[,c("Candidates","Gene2","Cancer_Type","Final_Score")])
    
    feature_df <- unique(reshape2::melt(feature_df,id.vars=c("Pair","Cancer_Type","Candidates")))
    score_df <- unique(reshape2::melt(score_df,id.vars=c("Pair","Cancer_Type","Candidates","Gene2")))
    final_score_df <- unique(reshape2::melt(final_score_df))
    colnames(feature_df)[4] <- "Feature"
    colnames(score_df)[5:6] <- c("Schema","Score")
    colnames(final_score_df)[4:5] <- c("Schema","Score")

    interesting_subset_feature <- feature_df %>% filter_at(vars(highlight_column), all_vars(. %in% highlight_list))
    interesting_subset_score <- score_df %>% filter_at(vars(highlight_column), all_vars(. %in% highlight_list))
    interesting_subset_final_score <- final_score_df %>% filter_at(vars(highlight_column), all_vars(. %in% highlight_list))
##    interesting_subset_score <- 

##    str(interesting_subset_feature)
##    str(score_df)
##    str(final_score_df)
    interesting_subset_final_score <- interesting_subset_final_score %>% group_by(Cancer_Type) %>% mutate(Rank=(length(Score)+1)-rank(Score)) %>% data.frame
    

    rownames(interesting_subset_final_score) <- interesting_subset_final_score$Candidates
    

    
    
    
    interesting_subset_feature$Rank <- interesting_subset_final_score[interesting_subset_feature$Candidates,"Rank"]
##    interesting_subset_feature <- interesting_subset_feature[unlist(lapply(unique(interesting_subset_feature$Pair), function(x) min(grep(x, interesting_subset_feature$Pair)))),]
    interesting_subset_feature$Pair <- paste0(interesting_subset_feature$Pair,": Rank ",interesting_subset_feature$Rank)

    interesting_subset_final_score$Gene2 <- paste0(interesting_subset_final_score$Gene2,": Rank ",interesting_subset_final_score$Rank) 
##    str(interesting_subset_score)
    

   
    p <- ggplot(score_df, aes(Score, fill=Cancer_Type))+geom_histogram(position="identity")+geom_vline(data=interesting_subset_score, aes(xintercept=Score, color=Cancer_Type),size=0.4,linetype="dashed")+facet_grid(Cancer_Type~Schema, scales="free")+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=12),strip.text.y = element_text(colour = "black", face = "bold",size=8))+scale_TF_score(scale,extension=".txt")+geom_text(data=interesting_subset_score,size=fontsize ,aes(x=Score-0.1,color=Cancer_Type ,y=+Inf,label=Pair,hjust=1, angle=90, fontface="bold"),check_overlap=FALSE)+guides(fill=guide_legend(ncol=1))


    p2_list <- lapply(unique(feature_df$Feature), function(x) dplyr::filter(feature_df, Feature==x))
    p2_list <- lapply(p2_list , function(x) ggplot(x, aes(value, fill=Cancer_Type))+geom_histogram(position="identity",aes(alpha=0.7))+geom_vline(data=dplyr::filter(interesting_subset_feature, Feature==unique(x$Feature)), aes(xintercept=value, color=Cancer_Type),size=0.4,linetype="dashed")+facet_wrap(Cancer_Type~.,scales="free",ncol=2)+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=20),axis.text.y=element_text(face='bold',size=20),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=12),strip.text.y = element_text(colour = "black", face = "bold",size=12),plot.title = element_text(size=30))+scale_TF_score(scale,extension=".txt")+geom_text_repel(data=dplyr::filter(interesting_subset_feature, Feature==unique(x$Feature)),size=2*fontsize ,aes(x=value,y=+Inf,label=Pair,hjust=1, angle=90, fontface="bold"),segment.size=0.25,check_overlap=FALSE)+guides(fill=guide_legend(ncol=1))+labs(title=unique(x$Feature)))


##    p2 <- ggplot(feature_df, aes(value, fill=Cancer_Type))+geom_histogram(position="identity")+geom_vline(data=interesting_subset_feature, aes(xintercept=value, color=Cancer_Type),size=0.4,linetype="dashed")+facet_grid(Cancer_Type~Feature,scales="free")+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=12),axis.text.y=element_text(face='bold',size=15),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=15),strip.text.y = element_text(colour = "black", face = "bold",size=8))+scale_TF_score(scale,extension=".txt")+geom_text_repel(data=interesting_subset_feature,size=fontsize ,aes(x=value,color=Cancer_Type ,y=+Inf,label=Pair,hjust=1, angle=90, fontface="bold"),segment.size=0.25,check_overlap=FALSE)+guides(fill=guide_legend(ncol=1))

    p3 <- ggplot(score_df, aes(Score, fill=Cancer_Type))+geom_histogram(position="identity")+geom_vline(data=interesting_subset_final_score, aes(xintercept=Score, color=Cancer_Type),size=0.4,linetype="dashed")+facet_wrap(Cancer_Type~., scales="free_x",nrow=2)+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=20),axis.text.y=element_text(face='bold',size=20),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=12),strip.text.y = element_text(colour = "black", face = "bold",size=8))+scale_TF_score(scale,extension=".txt")+geom_text_repel(data=interesting_subset_final_score,size=2*fontsize ,aes(x=Score-0.05,color=Cancer_Type ,y=+Inf,label=Gene2,hjust=1, angle=90, fontface="bold"),segment.size=0.25,check_overlap=FALSE)+guides(fill=guide_legend(ncol=1))


    print(p)
    lapply(p2_list, function(x) print(x))
##    print(p2)
    print(p3)

#return(interesting_subset_feature)
}




convert_long_predictions_to_matrix <- function(df, TSNE_df, normalized=TRUE, signif_only=TRUE){
    require(dplyr)
    if(signif_only==TRUE){ df <- dplyr::filter(df, Signif==TRUE); print("Only keeping significant predictions")} else{ df <- df}

    df <- df %>% group_by(Subtype,Prediction) %>% mutate(Freq=length(unique(Cell_Line))) %>% data.frame

    df2 <- df %>% group_by(Subtype) %>% mutate(Freq=round(Freq/length(unique(Cell_Line)),digits=2)) %>% data.frame

##    str(df)
##    str(df2)

    if(normalized==TRUE){
    out <- unique(df2[c("Subtype","Prediction","Freq")])} else if(normalized==FALSE){ out <- unique(df[c("Subtype","Prediction","Freq")])}

    mat <- matrix(0, nrow=length(unique(TSNE_df$Disease)), ncol=length(unique(TSNE_df$Subtype)), dimnames=list(unique(TSNE_df$Disease),unique(TSNE_df$Subtype)))

    print("Filling matrix")

##    str(out)

    for(i in 1:nrow(out)){ mat[out[i,2],out[i,1]] <- out[i,3]}
    return(mat)
}


parse_CMS_caller <- function(CMS_caller_output, qval_cutoff=0.05, Achilles_metadata){
    require(dplyr)
    colnames(CMS_caller_output)[1] <- "Prediction"
    CMS_caller_output$Signif <- CMS_caller_output$FDR <= qval_cutoff
    CMS_caller_output$Prediction <- as.character(CMS_caller_output$Prediction)
    CMS_caller_output$Cell_Line <- rownames(CMS_caller_output)
    CMS_caller_output$Subtype <- Achilles_metadata[rownames(CMS_caller_output),"TCGA"]
    CMS_caller_output$Detailed_Cancer <- Achilles_metadata[CMS_caller_output$Cell_Line,"lineage_subtype"]
    CMS_caller_output$Detailed_Cancer_Ext <- Achilles_metadata[CMS_caller_output$Cell_Line,"disease_subtype"]
    CMS_caller_output$Primary <- Achilles_metadata[CMS_caller_output$Cell_Line,"primary_or_metastasis"]
    CMS_caller_output <- CMS_caller_output %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
    
    CMS_caller_output <- CMS_caller_output%>% group_by(Subtype, Primary) %>% mutate(Accuracy_fraction=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
CMS_caller_output <- CMS_caller_output%>% group_by(Subtype) %>% mutate(Accuracy_fraction_overall=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
    CMS_caller_output <- CMS_caller_output%>% group_by(Subtype,Signif) %>% mutate(Accuracy_fraction_signif=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
    CMS_caller_output$Q_val_cutoff <- qval_cutoff
    CMS_caller_output$Num_Predictions <- 1

##    str(CMS_caller_output)


    cols <- ncol(CMS_caller_output)
  ##  print(cols)
  ##  print(colnames(CMS_caller_output))

    int <- unique(CMS_caller_output[c(1,26:cols)])

##    str(int)
    disease_types <- intersect(unique(unlist(lapply(int[,1], function(x) unlist(strsplit(x,"\\."))))),int$Subtype)
        
    ##str(disease_types)
    int <- dplyr::filter(int, Subtype %in% disease_types)

    return(int)

}
plot_CMS_caller <- function(parsed_CMS_caller_output,Signif_filter=TRUE, facet=NULL){
    require(dplyr)
    require(ggplot2)
    if(Signif_filter==TRUE){
        parsed_CMS_caller_output <- dplyr::filter(parsed_CMS_caller_output,Signif==TRUE)} else {parsed_CMS_caller_output <- parsed_CMS_caller_output}

  
    
    

    if(is.null(facet) == TRUE){

        int <- as.data.frame(table(parsed_CMS_caller_output[c("Prediction","Subtype")]),stringsAsFactors=F)
        int <- int %>% group_by(Subtype) %>% mutate(Freq_Norm=round(Freq/sum(Freq),digits=2)) %>% data.frame

    p <- ggplot(int, aes(y=Prediction,Subtype,fill=Freq,label=Freq))+geom_tile()+geom_text(fontface="bold",size=6,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Raw Numbers")
        p2 <- ggplot(int, aes(y=Prediction,Subtype,fill=Freq_Norm,label=Freq_Norm))+geom_tile()+geom_text(fontface="bold",size=4.5,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")
        p3 <- ggplot(out, aes(y=Prediction,Subtype,fill=Freq_Norm,label=Freq))+geom_tile()+geom_text(fontface="bold",size=4.5,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")
    } else {
        out <- data.frame(stringsAsFactors=F)
        for(i in unique(parsed_CMS_caller_output[,facet])) {
            print(i)
            sub <- parsed_CMS_caller_output[which(parsed_CMS_caller_output[,facet] ==i),]
            int <- as.data.frame(table(sub[c("Prediction","Subtype")]),stringsAsFactors=F)
            int <- int %>% group_by(Subtype) %>% mutate(Freq_Norm=round(Freq/sum(Freq),digits=2)) %>% data.frame
            int$Method <- i
            out <- rbind(out,int,stringsAsFactors=F)
        }                                                        

       
        p <- ggplot(out, aes(y=Prediction,Subtype,fill=Freq,label=Freq))+geom_tile()+geom_text(fontface="bold",size=6,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Raw Numbers")+facet_wrap(as.formula(paste0(facet,"~.")),ncol=2)
        p2 <- ggplot(out, aes(y=Prediction,Subtype,fill=Freq_Norm,label=Freq_Norm))+geom_tile()+geom_text(fontface="bold",size=4.5,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")+facet_wrap(as.formula(paste0(facet,"~.")),ncol=2)

        p3 <- ggplot(out, aes(y=Prediction,Subtype,fill=Freq_Norm,label=Freq))+geom_tile()+geom_text(fontface="bold",size=4.5,color="gray50")+scale_fill_viridis_c(option="magma")+theme(axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+labs(title="Signif Predictions: Normalized Numbers")+facet_wrap(as.formula(paste0(facet,"~.")),ncol=2)}
    

    print(p)
    print(p2)
    print(p3)
    }
    
plot_key_TF_essentiality <- function(Achilles_df, NTP_results,TSNE_df,driver_column=c("Cluster_drivers","Signif_drivers","Top_5")){
    require(reshape2)
    require(dplyr)
    require(ggsignif)

    rownames(NTP_results) <- NTP_results[,1]
    TF_df <- unique(TSNE_df[c("Disease",driver_column)])
    rownames(TF_df) <- TF_df[,1]
##    str(TF_df)

    Achilles_df <- Achilles_df[,intersect(colnames(Achilles_df),NTP_results[,1])]
    Achilles_df <- reshape2::melt(as.matrix(Achilles_df))

    Achilles_df[1:2] <- apply(Achilles_df[1:2],2,as.character)
    colnames(Achilles_df)[1:2] <- c("Gene","Cell_Line")
    Achilles_df$Subtype <- NTP_results[Achilles_df[,2],"Subtype"]
    Achilles_df$Prediction <- NTP_results[Achilles_df[,2],"Prediction"]
    Achilles_df$Signif <- NTP_results[Achilles_df[,2],"Signif"]
    Achilles_df <- Achilles_df %>% group_by(Prediction,Gene) %>% mutate(search_term=paste0("^",Gene,"-|-",Gene,"-|-",Gene,"$")) %>%  mutate(Key_TF= grepl(search_term, TF_df[Prediction,2])) %>% data.frame
##    str(Achilles_df)
    ##    str(filter(Achilles_df, Key_TF==TRUE))
    pan_can_df <- Achilles_df %>% group_by(Gene) %>% dplyr::select(value) %>% summarize_all(median) %>% ungroup
  ##  str(pan_can_df)
  ##  print(table(pan_can_df$value <=-0.5))
    pan_can_genes <- dplyr::filter(pan_can_df,value <= -0.5)$Gene
  ##  str(pan_can_genes)
    Achilles_df$Pan_Can_essential <- Achilles_df$Gene %in% pan_can_genes
    print(table(Achilles_df[c("Pan_Can_essential","Key_TF")]))

    Achilles_df$Flag <- paste0(Achilles_df$Key_TF,"~",Achilles_df$Pan_Can_essential)
  ##  print(table(Achilles_df$Flag))
    flag_translator <- list("TRUE~TRUE"="Key_TF","TRUE~FALSE"="Key_TF","FALSE~TRUE"="Pan_can_essential","FALSE~FALSE"="Other")
    Achilles_df$Key_TF <- unlist(flag_translator[Achilles_df$Flag])
  ##  print(table(Achilles_df$Key_TF))


    sub_Achilles_df <- Achilles_df %>% group_by(Prediction,Gene,Key_TF) %>% dplyr::select(value) %>% summarize_all(mean) %>% ungroup

##    str(sub_Achilles_df)
##    print(table(sub_Achilles_df[c("Key_TF","Prediction")]))


    mean_df <- Achilles_df %>% group_by(Prediction,Key_TF) %>% dplyr::select(value) %>% summarize_all(median) %>% ungroup
    ##    str(mean_df)

    p <- ggplot(Achilles_df,aes(value, fill=Key_TF,alpha=0.3))+geom_density(adjust=0.5)+facet_wrap(Prediction~.,ncol=2)+labs(title="Gene essentiality for all cell lines, critical TFs highlighted")+theme(axis.text.x=element_text(face='bold',size=10,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=6))+geom_vline(xintercept=-0.5)+geom_vline(data=mean_df,aes(alpha=0.3,xintercept=value,color=Key_TF),linetype=2, size=1)+scale_fill_manual(values=c("red","gray50","cyan"))+scale_color_manual(values=c("red","gray50","cyan"))
    p2 <- ggplot(Achilles_df,aes(value, fill=Key_TF,alpha=0.3))+geom_density(adjust=0.5)+facet_wrap(Subtype~.,ncol=2)+labs(title="Gene essentiality for all cell lines, critical TFs highlighted")+theme(axis.text.x=element_text(face='bold',size=10,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=6))+geom_vline(xintercept=-0.5)+scale_fill_manual(values=c("red","gray50","cyan"))

    p3 <- ggplot(sub_Achilles_df,aes(Key_TF, value, fill=Key_TF,alpha=0.3))+geom_violin()+facet_wrap(Prediction~.,ncol=ncol)+labs(title="Gene essentiality for all cell lines, critical TFs highlighted")+theme(axis.text.x=element_text(face='bold',size=10,hjust=0.95,angle=45),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=6))+geom_hline(yintercept=-0.5)+scale_fill_manual(values=c("red","gray50","cyan"))+geom_signif(comparisonsb=list(c("Key_TF","Other")), y_position = max(sub_Achilles_df$value)-0.1,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25), textsize=5,test="t.test")
    

    print(p)
    print(p2)
    print(p3)
    ##str(filter(Achilles_df, Key_TF==TRUE))
    ##print(table(Achilles_df$Key_TF))
}

rowMedians <- function(df,na.rm=TRUE){

    out <- apply(df,1,function(x) median(x, na.rm=na.rm))
    names(out) <- rownames(df)
    return(out)}

make_fisher_df_key_TFs <- function(top50_list, TSNE_df, num_TFs=10,Achilles_df,group_column="Disease",TF_column="Signif_drivers"){
    TF_df <- unique(TSNE_df[c(group_column,TF_column)])
    TF_df_list <- lapply(TF_df[,2], function(x) unlist(strsplit(x,"-")))
    names(TF_df_list) <- TF_df[,1]
    
    key_TFs_in_top50 <- unlist(lapply(names(top50_list), function(x) length(intersect(top50_list[[x]],TF_df_list[[x]][1:num_TFs]))))

    fisher_df <- as.data.frame(cbind(key_TFs_in_top50,num_TFs, names(top50_list)),stringsAsFactors=F)
    fisher_df[1:2] <- apply(fisher_df[1:2],2,function(x) as.numeric(x))
    colnames(fisher_df) <- c("Num_Key_TFs","TF_sample_size","Cluster")

##    str(fisher_df)
    vec <- lapply(1:nrow(fisher_df), function(x) c(rep("Yes",fisher_df[x,1]),rep("No",fisher_df[x,2]-fisher_df[x,1])))
    names(vec) <- fisher_df[,3]
##    str(vec)
    universe <- lapply(names(top50_list), function(x) c(rep("Yes",length(top50_list[[x]])),rep("No",nrow(Achilles_df)-length(top50_list[[x]]))))
    names(universe) <- names(top50_list)

    out <- list(vec,universe,fisher_df)
    return(out)

}
    
TF_score_vs_essentiality <- function(achilles_df,NTP_df ,TSNE_df,TF_list="Genome_Wide_TFs.txt",Essential_list="Depmap_essential_genes.txt"){

    require(reshape2)
    require(dplyr)
    achilles_df <- achilles_df[intersect(colnames(achilles_df),NTP_df$Cell_Line)]
    melted_df <- reshape2::melt(as.matrix(achilles_df))
    colnames(melted_df) <- c("Gene","Cell_Line","Essentiality")
    
    melted_df$Prediction <- NTP_df[as.character(melted_df$Cell_Line),"Prediction"]
    ##melted_df[1:2]  <- apply(melted_df[1:2],2,function(x) as.character(x))

    TF_df <- unique(TSNE_df[c("Disease","Signif_drivers")])
    rownames(TF_df) <- TF_df[,1]
    TF_df_list <- lapply(TF_df[,2], function(x) unlist(strsplit(x,"-")))
    names(TF_df_list) <- TF_df[,1]
    ###str(TF_df_list)
    ##str(melted_df)

    melted_df <- melted_df %>% group_by(Prediction,Gene) %>% mutate(search_term=paste0("^",Gene,"-|-",Gene,"-|-",Gene,"$")) %>% mutate(Key_TF=grepl(unique(search_term), TF_df[Prediction,2])) %>% data.frame
    
  ##  print(table(melted_df$Key_TF))
    ##str(melted_df)
    index <- which(melted_df$Key_TF==TRUE)
    melted_df$TF_score_rank <- 1000
    if(length(index) >=1){
    melted_df[index,"TF_score_rank"] <- unlist(lapply(index, function(x) grep(paste0(melted_df[x,"Gene"],"$"), TF_df_list[[melted_df[x,"Prediction"]]])))} else{ print("None of the key TFs are in this dataset")}
    ##str(melted_df)

    ##print(table(melted_df$TF_score_rank))
    summary_df <- melted_df %>% group_by(Prediction,Gene) %>% mutate(Mean_Essentiality=mean(Essentiality,na.rm=TRUE), Median_Essentiality=median(Essentiality,na.rm=TRUE),SD_Essentiality=sd(Essentiality,na.rm=TRUE)) %>% mutate(Median_Upper=Median_Essentiality+SD_Essentiality, Median_Lower=Median_Essentiality-SD_Essentiality,Mean_Upper=Mean_Essentiality+SD_Essentiality,Mean_Lower=Mean_Essentiality-SD_Essentiality) %>% data.frame
    summary_df <- unique(summary_df[c(-2,-3,-5)])

    if(file.exists(TF_list) ==TRUE){ TF_list <- readLines(TF_list)
                                     melted_df$Is_TF <- melted_df$Gene %in% TF_list
                                     summary_df$Is_TF <- summary_df$Gene %in% TF_list

                                 } else{print("We don't have a list of TFs in this dataset, not flagging")}

    if(file.exists(Essential_list) ==TRUE){ Essential_list <- readLines(Essential_list)
                                     melted_df$Pan_can_essential <- melted_df$Gene %in% Essential_list
                                     summary_df$Pan_can_essential <- summary_df$Gene %in% Essential_list

                                 } else{print("We don't have a list of TFs in this dataset, not flagging")}

##    str(summary_df)
    return(list(melted_df,summary_df))
}


plot_TF_score_vs_essentiality <- function(TF_score_essentiality_output, rank_cutoff=1000,scale="/pbtech_mounts/homes024/anf2034/ANF_Final_Cluster_key.txt"){
    require(ggplot2)
    require(dplyr)
    require(ggrepel)

    summary_df <- TF_score_essentiality_output[[2]]
    if(length(grep("Pan_can",colnames(summary_df)))>0){

    p <- ggplot(dplyr::filter(summary_df,TF_score_rank<rank_cutoff), aes(TF_score_rank,Mean_Essentiality,color=Prediction,label=Gene,shape=Pan_can_essential))+geom_point(size=4)+geom_errorbar(aes(ymin=Mean_Lower, ymax=Mean_Upper))+geom_smooth(method="lm")+facet_wrap(scales="free",Prediction~., ncol=4)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_hline(yintercept=-0.5,size=2, alpha=0.3, linetype=2)+scale_y_reverse()+scale_x_reverse()+geom_text_repel(size=2,fontface="bold", segment.size=0.25)+labs(title="Mean Essentiality")+scale_TF_score(scale)
    p2  <- ggplot(dplyr::filter(summary_df,TF_score_rank<rank_cutoff), aes(TF_score_rank,Median_Essentiality,color=Prediction,label=Gene,shape=Pan_can_essential))+geom_point(size=4)+geom_errorbar(aes(ymin=Median_Lower, ymax=Median_Upper))+geom_smooth(method="lm")+facet_wrap(scales="free",Prediction~., ncol=4)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_hline(yintercept=-0.5,size=2, alpha=0.3, linetype=2)+scale_y_reverse()+scale_x_reverse()+geom_text_repel(size=2,fontface="bold", segment.size=0.25)+labs(title="Median Essentiality")+scale_TF_score(scale)
} else{
            p <- ggplot(dplyr::filter(summary_df,TF_score_rank<rank_cutoff), aes(TF_score_rank,Mean_Essentiality,color=Prediction,label=Gene))+geom_point(size=4)+geom_errorbar(aes(ymin=Mean_Lower, ymax=Mean_Upper))+geom_smooth(method="lm")+facet_wrap(scales="free",Prediction~., ncol=4)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_hline(yintercept=-0.5,size=2, alpha=0.3, linetype=2)+scale_y_reverse()+scale_x_reverse()+geom_text_repel(size=2,fontface="bold", segment.size=0.25)+labs(title="Mean Essentiality")+scale_TF_score(scale)
    p2  <- ggplot(dplyr::filter(summary_df,TF_score_rank<rank_cutoff), aes(TF_score_rank,Median_Essentiality,color=Prediction,label=Gene))+geom_point(size=4)+geom_errorbar(aes(ymin=Median_Lower, ymax=Median_Upper))+geom_smooth(method="lm")+facet_wrap(scales="free",Prediction~., ncol=4)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_hline(yintercept=-0.5,size=2, alpha=0.3, linetype=2)+scale_y_reverse()+scale_x_reverse()+geom_text_repel(size=2,fontface="bold", segment.size=0.25)+labs(title="Median Essentiality")+scale_TF_score(scale)
        }
   

    print(p)
    print(p2)

}


FPKM_to_TPM <- function(fpkm) { out <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
                                return(out)}

FPKM_to_TPM_mat <- function(mat) { mat <- as.data.frame(apply(mat,2, function(x) log2(FPKM_to_TPM(x))))

                                   return(mat)}
DE_gene_from_TPM <- function(mat, metadata, reference_level,column,label="Sample"){
    require(dplyr)
    require(limma)

### metadata df requires three columns minimum, sample ID, condition (case control etc) and the column name in the original data the ID corresponds to (Label)  [ID,Condition,Label]. Additional metadata columns can be present but are ignored.
    
    metadata[,column] <- relevel(factor(metadata[,column]), reference_level)
    

    print("Checking to see if metadata and data match up")
    if(ncol(mat) == nrow(metadata)) { print("Same number of conditions detected, please double-check that they're ordered in the same way")} else{ print("The metadata and data don't match, subsetting data")

                                                                                                                                               mat2 <- mat[,intersect(metadata[,label],colnames(mat))]
                                                                                                                                               print(paste0("Keeping ", length(intersect(metadata[,label],colnames(mat)))," out of ",ncol(mat)," columns from original matrix"))}
    metadata <- metadata[intersect(metadata[,label],colnames(mat2)),]
    design <- model.matrix(~metadata[,column])
    fit <- lmFit(mat2,design)
    comparisons <- paste0(levels(metadata[,column])[-1],"_vs_",levels(metadata[,column])[1])
    coeffs <- 2:nlevels(metadata[,column])

    ##print(comparisons)
    ##print(coeffs)

    out <- topTable(eBayes(fit),number=nrow(mat))
    colnames(out)[1:length(coeffs)] <- paste0(comparisons,"_FC")
##    print(head(out))
    out$Gene <- rownames(out)

    pval_list <- lapply(coeffs, function(x) topTable(eBayes(fit),coef=x,number=nrow(mat))[c("t","P.Value","adj.P.Val")])
##    str(pval_list)
##    final_pval_df <- data.frame(stringsAsFactors=F)

    for(i in 1:length(pval_list)){ sub <- pval_list[[i]]
                                   colnames(sub) <- paste0(comparisons[i],"_",colnames(sub))
                                   sub$Gene <- rownames(sub)
                                   pval_list[[i]] <- sub
                               }
    final_pval_df <- Reduce(function(...) merge(...,all=TRUE, by="Gene"),pval_list)
##    str(final_pval_df)
    out <- merge(out,final_pval_df, by="Gene")
    ##str(out)
    rownames(out) <- out$Gene
    return(out)
    
    
}
   

 
DGTAC_validation <- function(TF_bed,peak_gene_bed,TF, TF_targets,TF_cutoff=700,enhancer_gr="DGTAC_validation/distal_regions.bed"){
    require(dplyr)
    require(sfsmisc)
    

    TF_bed <- bed_to_granges_dynamic(TF_bed,header=FALSE)
    TF_bed <- TF_bed[which(TF_bed$V5 >=TF_cutoff)]
    
    ##print(TF_bed)
    peak_gene_bed <- bed_to_granges_dynamic(peak_gene_bed,header=FALSE)
    ##print(peak_gene_bed)
    if(is.null(enhancer_gr)==FALSE){
        enhancer_gr <- bed_to_granges_dynamic(enhancer_gr)
        peak_gene_bed <- intersect_with_metadata(peak_gene_bed, enhancer_gr)} else{ peak_gene_bed <- peak_gene_bed}
    ##print(peak_gene_bed)

    common <- intersect(peak_gene_bed$V5, TF_targets)

    total <- union(peak_gene_bed$V5, TF_targets)

    ##str(common)
    ##str(total)

    vec <- seq(0,1,0.05)

    out <- data.frame(stringsAsFactors=F)
    for(i in vec){
        sub <- peak_gene_bed[which(peak_gene_bed$V6 >=i)]
        ##print(length(sub))
        predicted_targets <- intersect_with_metadata(sub,TF_bed)
        predicted_targets <- unique(predicted_targets$V5)
        ##str(intersect(predicted_targets,TF_targets))
        TP <- length(intersect(predicted_targets,TF_targets))
        FP <- length(setdiff(predicted_targets,TF_targets))
        FN <- length(setdiff(predicted_targets,x=TF_targets))
        if(TP>0){
            Precision <- TP/(TP+FP)
            Recall <- TP/(TP+FN)} else { Precision <- 1
                                         Recall <-0}
        out <- rbind(out, cbind(TP,FP,FN,Precision,Recall,i,TF))
    }

    out[1:6] <- apply(out[1:6],2,as.numeric)
    colnames(out)[6] <- "Cutoff"
    out[is.na(out)] <- 1
    if(all(out$TP==0)==TRUE){ out$Area <- 0} else{    out$Area <- integrate.xy(out$Recall, out$Precision)}

    out$Random_expectation <- (out$TP+out$FN)/sum(out[1,1:3])
    
    return(out)}

translate_ENSMBL_to_HGNC <- function(ENSEMBL_IDs,version="EnsDb.Hsapiens.v79",key="GENEID", column="GENENAME",split=TRUE){
    require(version,character.only=T)
    if(split==TRUE){ENSEMBL_IDs <- unlist(lapply(ENSEMBL_IDs, function(x) unlist(strsplit(x,"\\."))[1]))} else { ENSEMBL_IDs <- ENSEMBL_IDs}

    ENSEMBL_IDs <- unique(ENSEMBL_IDs)
##    str(ENSEMBL_IDs)
    symbols <- ensembldb::select(get(version),keys=ENSEMBL_IDs,keytype=key, columns=column)
    rownames(symbols)  <- symbols[,1]
    final_symbols <- symbols[ENSEMBL_IDs,2]
    return(final_symbols)}

make_and_run_fishers_exact_vectors <- function(vec, universe, truth_set,label=NULL,alternative_opt=c("two.sided","greater","less"),log_OR=TRUE){
    vec_int <- vec %in% truth_set
    universe_int <- universe %in% truth_set

##    str(vec_int)
##    str(universe_int)

##    print(table(vec_int))
##    print(table(universe_int))
    out <- fishers_exact_vec(vec_int,universe_int,label, alternative_opt=alternative_opt, log_OR=log_OR)

    return(out)}

plot_fisher_output <- function(fisher_output,title=NULL,shape=NULL,size=5){
    require(ggplot2)
    x_lab <- ifelse(unique(fisher_output$Is_Log_Odds)==TRUE,"Log Odds Ratio","Odds Ratio")
    hline <- ifelse(unique(fisher_output$Is_Log_Odds)==TRUE,0,1)
    fisher_output$Alpha <- 1-fisher_output$P_value
if(is.null(shape)){
    p <- ggplot(fisher_output, aes(Odds_Ratio,Label, fill=Label))+geom_point(size=size,shape=21,color="gray50",stroke=0.1)+xlab(x_lab)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_errorbarh(aes(xmin=OR_Lower,xmax=OR_Upper,color=Label),height=0.05,size=1)+labs(title=title)+geom_vline(xintercept=hline, size=1, alpha=0.2, linetype=2)} else{

        p <- ggplot(fisher_output, aes_string("Odds_Ratio","Method", color="Label",shape=shape))+geom_point(size=size)+xlab(x_lab)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text.x = element_text(colour = "black", face = "bold",size=10),strip.text.y = element_text(colour = "black", face = "bold",size=8,angle=360))+geom_errorbarh(aes(xmin=OR_Lower,xmax=OR_Upper),height=0.05,size=1)+labs(title=title)+geom_vline(xintercept=hline, size=1, alpha=0.2, linetype=2)}

    return(p)}

plot_superimposed_graphs <- function(graph1,graph2,label1=NULL,label2=NULL,vsize=3,vcex=0.1){
    require(igraph)
    require(dplyr)
    combo_graph <- igraph::union(graph1,graph2,byname=TRUE)
    layg_final <- layout_in_circle(combo_graph)
                                         
    layg1 <- layout_in_circle(graph1)
    layg2 <- layout_in_circle(graph2)

#    layg2[which(V(graph2)$name %in% V(graph1)$name), ] <- layg1[which(V(graph1)$name %in% V(graph2)$name),]

    layg2[which(V(graph2)$name %in% V(combo_graph)$name), ] <- layg_final[which(V(combo_graph)$name %in% V(graph2)$name),]
    layg1[which(V(graph1)$name %in% V(combo_graph)$name), ] <- layg_final[which(V(combo_graph)$name %in% V(graph1)$name),]

    el1 <- as.data.frame(get.edgelist(graph1))
    el1$E <- paste0(el1[,1],"_",el1[,2])

    el2 <- as.data.frame(get.edgelist(graph2))
    el2$E <- paste0(el2[,1],"_",el2[,2])

    common1 <- which(el1$E %in% el2$E)
    common2 <- which(el2$E %in% el1$E) 

    

    E(graph2)$color <- adjustcolor("yellow", alpha.f = 0.3)
    E(graph1)$color <- adjustcolor("cyan",alpha.f=0.3)

    
   ## plot(graph1 , vertex.size=2, layout=layg1, rescale=FALSE,vertex.label.cex=0.05,width=2,edge.arrow.size=0.2)
   ## plot(graph2 , vertex.size=2, layout=layg2, rescale=FALSE,vertex.label.cex=0.05,width=2,edge.arrow.size=0.2)


    E(graph2)$color[common2] <- adjustcolor("purple", alpha.f = 0.3)
    E(graph1)$color[common1] <- adjustcolor("purple", alpha.f = 0.3)
    
    plot(graph1 , vertex.size=vsize, layout=layg1, rescale=FALSE,vertex.label.cex=vcex,vertex.label.color="black",width=0.05,edge.arrow.size=0.1,vertex.color="white")
    plot(graph2 , vertex.size=vsize, layout=layg2, rescale=FALSE,vertex.label.cex=vcex,vertex.label.color="black",width=0.05,add=T,edge.arrow.size=0.1,vertex.color="white")
    if(is.null(c(label1,label2))==TRUE){
        legend("bottomleft", inset=0.1, title="Networks",c("Graph1 only","Graph2 only","Common"), fill=c("cyan","yellow","purple"), horiz=TRUE, pt.cex=50,cex=50)} else if(is.null(c(label1,label2)) ==FALSE){
            legend("bottomleft", inset=0, title="Networks",c(paste0(label1," only"),paste0(label2," only"),"Common"), fill=c("cyan","yellow","purple"), horiz=TRUE, cex=0.5)} 
        ##layout=c("circle","circular","drl","fr","fr3d","fr_3d","grid_fr","kk","kk3d","kk_3d","large","lgl","large_graph","random","random_3d","rt","tree","rt_circular","sphere","spherical","circular_3d")   
}

get_quantile <- function(vec,abs_val=FALSE){
    if(abs_val==TRUE){
        out <- ecdf(abs(vec))(abs(vec))} else{ out <- ecdf(vec)(vec)}
    if(is.null(names(vec))== FALSE){ names(out) <- names(vec)} else{ print("")}
    
    return(out)}

DGTAC_validation_Fishers <- function(TF_bed,peak_gene_bed,TF, truth_set,label=NULL,alternative_opt=c("two.sided","greater","less"),log_OR=TRUE,TF_cutoff=700,link_cutoff=0.5,enhancer_gr="DGTAC_validation/distal_regions.bed",universe=NULL){
    
    require(dplyr)
    require(sfsmisc)
    

    TF_bed <- bed_to_granges_dynamic(TF_bed,header=FALSE)
    TF_bed <- TF_bed[which(TF_bed$V5 >=TF_cutoff)]
        
    peak_gene_bed <- bed_to_granges_dynamic(peak_gene_bed,header=FALSE)
    peak_gene_bed <- peak_gene_bed[which(abs(peak_gene_bed$V6) >=link_cutoff)]
    
    if(is.null(enhancer_gr)==FALSE){
        enhancer_gr <- bed_to_granges_dynamic(enhancer_gr)
        peak_gene_bed <- intersect_with_metadata(peak_gene_bed, enhancer_gr)} else{ peak_gene_bed <- peak_gene_bed}
    ##print(peak_gene_bed)

    predicted <- unique(intersect_with_metadata(peak_gene_bed,TF_bed)$V5)
    if(is.null(universe)==TRUE){universe <- unique(peak_gene_bed$V5)} else{
        ##str(universe);universe <- union(universe,peak_gene_bed$V5);str(universe)
        universe <- universe
        
                                                                        }
    truth_set <- intersect(truth_set, universe)
##    str(predicted)
##    predicted <- intersect(predicted,universe)

##    print(TF_bed)
##    print(peak_gene_bed)
##    str(predicted)
##    str(universe)
##    str(truth_set)

    out <- make_and_run_fishers_exact_vectors(predicted,universe,truth_set,label,alternative_opt=alternative_opt)
    if(nrow(out)>1){    out <- dplyr::filter(out, Category==TRUE);    out$Evaluation <- alternative_opt} else{ out <- out}
    return(out)}





PIQ_motifmatchr_comparison <- function(PIQ_output_dir,TF,motif_match_mat,TF_bed, peak_regions,peak_filter=-0.5,buffer=0,shiftb=0){
    require(sfsmisc)
    require(pracma)
    require(GenomicRanges)
    suffix <- c("DN","D","I","IN")
    search_terms <- paste0(paste0(TF,suffix),collapse="[0-9]?|")

    file_list <- list.files(PIQ_output_dir,".bed")

    TF_bedlist <- grep(search_terms, file_list,value=T)
    print(TF_bedlist)

    print(paste0("Found ", length(TF_bedlist)/2," bed files for ",TF))

    final_bed_list <- lapply(TF_bedlist, function(x) PIQ_bed_to_gr(paste0(PIQ_output_dir,x)))
    

    print("Finished importing bed files")


    peak_regions <- peak_regions[which(peak_regions$score >=peak_filter)]
##    print(peak_regions)
    PIQ_bed <- sort_gr(do.call("c",final_bed_list))
    PIQ_bed <- PIQ_bed+buffer
##    print(PIQ_bed)
    PIQ_bed <- shift(intersect_with_metadata(PIQ_bed,peak_regions),shiftb)
##    print(TF_bed)
    TF_bed <- unique(intersect_with_metadata(peak_regions,TF_bed))
    
    print("Finished subsetting")



    
    PIQ_sequence <- sort(c(700,seq(min(PIQ_bed$purity),1000,length.out=19)))
    Peak_sequence <- seq(-0.5, max(peak_regions$score),length.out=20)

    ##str(PIQ_sequence)
    ##str(Peak_sequence)


    PIQ_out <- data.frame(stringsAsFactors=F)
    PIQ_out <- rbind(PIQ_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(PIQ_out) <- c("TP","FP","FN","TN","i")

    for(i in PIQ_sequence){
        sub <- PIQ_bed[which(PIQ_bed$purity >=i)]
        sub <- intersect_with_metadata(peak_regions,sub)
        int <- intersect_with_metadata(sub,TF_bed)
        #print(int)
        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        PIQ_out <- rbind(PIQ_out,cbind(TP,FP,FN,TN,i))
    }
    colnames(PIQ_out) <- c("TP","FP","FN","TN","Cutoff")
    index <- which(PIQ_out$TP == 0)
    PIQ_out$Precision <- PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FP)
    PIQ_out$Precision[index] <- 0
    PIQ_out$Recall <-  PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FN)
        
    PIQ_out$Baseline <- length(TF_bed)/length(peak_regions)
    PIQ_out$Area <- abs(trapz(PIQ_out$Recall, PIQ_out$Precision))
    PIQ_out$TF <- TF
    PIQ_out$Method <- "PIQ"

  ##  str(PIQ_out)

    print("Finished PIQ")

    motif_match_out <- data.frame(stringsAsFactors=F)
    motif_match_out <- rbind(motif_match_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(motif_match_out) <- c("TP","FP","FN","TN","i")
    
    for(i in Peak_sequence){
        sub <- peak_regions[which(peak_regions$score >=i)]
        sub$Motif <- motif_match_mat[names(sub),TF]
        sub2 <- sub[which(sub$Motif==1)]
        
        int <- unique(intersect_with_metadata(sub2,TF_bed))

        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub2),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub2),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        motif_match_out <- rbind(motif_match_out, cbind(TP,FP,FN,TN,i))
    }
    colnames(motif_match_out) <- c("TP","FP","FN","TN","Cutoff")

    index <- which(motif_match_out$TP == 0)
    motif_match_out$Precision <- motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FP)
    motif_match_out$Precision[index] <- 0
    motif_match_out$Recall <-  motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FN)
    
        
    motif_match_out$Baseline <- length(TF_bed)/length(peak_regions)
    motif_match_out$Area <- abs(trapz(motif_match_out$Recall, motif_match_out$Precision))
    motif_match_out$TF <- TF
    motif_match_out$Method <- "Motifmatchr"

##    str(motif_match_out)
    print("Finished motif match")
    print(paste0("MotifMatch Area:",unique(motif_match_out$Area)))
    print(paste0("PIQ Area:",unique(PIQ_out$Area)))

    PIQ_out$index <- grep(700,PIQ_out$Cutoff)
    motif_match_out$index <- grep(-0.5, motif_match_out$Cutoff)+nrow(PIQ_out)
     

    
final_out <- rbind(PIQ_out,motif_match_out)
return(final_out)
}

PIQ_motifmatchr_comparison2 <- function(PIQ_output_dir,TF,motif_match_mat,TF_bed, peak_regions,peak_filter=-0.5,buffer=0,shiftb=0,purity_cutoff=700,PIQ_bed=NULL){
    require(sfsmisc)
    require(pracma)
    require(GenomicRanges)

    ##    print(PIQ_bed)
    peak_regions <- peak_regions[which(peak_regions$score >=peak_filter)]


    if(is.null(PIQ_bed) == TRUE){
    suffix <- c("DN","D","I","IN")
    search_terms <- paste0(paste0(TF,suffix),collapse="[0-9]?|")

    file_list <- list.files(PIQ_output_dir,".bed")

    TF_bedlist <- grep(search_terms, file_list,value=T)
    print(TF_bedlist)

    print(paste0("Found ", length(TF_bedlist)/2," bed files for ",TF))

    final_bed_list <- lapply(TF_bedlist, function(x) PIQ_bed_to_gr(paste0(PIQ_output_dir,x)))


    print("Finished importing bed files")

##    print(length(peak_regions))
##    peak_regions <- peak_regions[which(peak_regions$score >=peak_filter)]
##    print(length(peak_regions))
    PIQ_bed <- sort_gr(do.call("c",final_bed_list))} else { PIQ_bed <- PIQ_bed}
    
    PIQ_bed <- PIQ_bed+buffer
    PIQ_bed <- GenomicRanges::shift(PIQ_bed,shiftb)

##    PIQ_bed <- shift(intersect_with_metadata(PIQ_bed,peak_regions),shiftb)
    PIQ_bed <- PIQ_bed[which(PIQ_bed$purity >=purity_cutoff)]
##    print(PIQ_bed)

    ##    print(TF_bed)

    if(is.character(TF_bed) ==TRUE){ TF_bed <- bed_to_granges(TF_bed)} else { TF_bed <- TF_bed}
    TF_bed <- unique(intersect_with_metadata(peak_regions,TF_bed))
    
    print("Finished subsetting")



    
   ## PIQ_sequence <- sort(c(700,seq(min(PIQ_bed$purity),1000,length.out=19)))
    Peak_sequence <- sort(c(-0.5,seq(min(peak_regions$score), max(peak_regions$score),length.out=20)))

    ##str(PIQ_sequence)
    ##str(Peak_sequence)


    PIQ_out <- data.frame(stringsAsFactors=F)
    PIQ_out <- rbind(PIQ_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(PIQ_out) <- c("TP","FP","FN","TN","i")

    for(i in Peak_sequence){
        sub <- peak_regions[which(peak_regions$score >=i)]
        sub <- intersect_with_metadata(sub,PIQ_bed)
        int <- intersect_with_metadata(sub,TF_bed)
        #print(int)
        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        PIQ_out <- rbind(PIQ_out,cbind(TP,FP,FN,TN,i))
    }
    colnames(PIQ_out) <- c("TP","FP","FN","TN","Cutoff")
    index <- which(PIQ_out$TP == 0)
    PIQ_out$Precision <- PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FP)
    PIQ_out$Precision[index] <- 0
    PIQ_out$Recall <-  PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FN)
        
    PIQ_out$Baseline <- length(TF_bed)/length(peak_regions)
    PIQ_out$Area <- abs(trapz(PIQ_out$Recall, PIQ_out$Precision))
    PIQ_out$TF <- TF
    PIQ_out$Method <- "PIQ"

  ##  str(PIQ_out)

    print("Finished PIQ")

    motif_match_out <- data.frame(stringsAsFactors=F)
    motif_match_out <- rbind(motif_match_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(motif_match_out) <- c("TP","FP","FN","TN","i")
    
    for(i in Peak_sequence){
        sub <- peak_regions[which(peak_regions$score >=i)]
        sub$Motif <- motif_match_mat[names(sub),TF]
        sub2 <- sub[which(sub$Motif==1)]
        
        int <- unique(intersect_with_metadata(sub2,TF_bed))

        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub2),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub2),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        motif_match_out <- rbind(motif_match_out, cbind(TP,FP,FN,TN,i))
    }
    colnames(motif_match_out) <- c("TP","FP","FN","TN","Cutoff")

    index <- which(motif_match_out$TP == 0)
    motif_match_out$Precision <- motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FP)
    motif_match_out$Precision[index] <- 0
    motif_match_out$Recall <-  motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FN)
    
        
    motif_match_out$Baseline <- length(TF_bed)/length(peak_regions)
    motif_match_out$Area <- abs(trapz(motif_match_out$Recall, motif_match_out$Precision))
    motif_match_out$TF <- TF
    motif_match_out$Method <- "Motifmatchr"

##    str(motif_match_out)
    print("Finished motif match")
    print(paste0("MotifMatch Area for ",TF,":",round(unique(motif_match_out$Area),2)))
    print(paste0("PIQ Area for ",TF,":",round(unique(PIQ_out$Area),2)))

    PIQ_out$index <- grep(-0.5,PIQ_out$Cutoff)
    motif_match_out$index <- grep(-0.5, motif_match_out$Cutoff)+nrow(PIQ_out)
     
    final_out <- rbind(PIQ_out,motif_match_out)
    final_out$Difference <- final_out$Area-final_out$Baseline
return(final_out)
}

HINT_motifmatchr_comparison <- function(HINT_bed,TF,motif_match_mat,TF_bed, peak_regions,peak_filter=-0.5,buffer=0,shiftb=0,purity_cutoff=5){
    require(sfsmisc)
    require(pracma)
    require(GenomicRanges)

    if(is.character(HINT_bed)==TRUE){
        print("working 1")
        PIQ_bed <- bed_to_granges_dynamic(HINT_bed)} else if(class(HINT_bed) == "GRanges"){
            print("working 2")
            PIQ_bed <- HINT_bed}
    else{print("Nothing")}
    PIQ_bed <- PIQ_bed[which(PIQ_bed$V5 >=purity_cutoff),]
##    print(PIQ_bed)

    print("Finished importing bed files")


    peak_regions <- peak_regions[which(peak_regions$score >=peak_filter)]
##    print(peak_regions)
    TF_bed <- bed_to_granges(TF_bed)
    TF_bed <- unique(intersect_with_metadata(peak_regions,TF_bed))
    
    print("Finished subsetting")



    
   ## PIQ_sequence <- sort(c(700,seq(min(PIQ_bed$purity),1000,length.out=19)))
    Peak_sequence <- sort(c(-0.5,seq(min(peak_regions$score), max(peak_regions$score),length.out=15)))

    ##str(PIQ_sequence)
    ##str(Peak_sequence)


    PIQ_out <- data.frame(stringsAsFactors=F)
    PIQ_out <- rbind(PIQ_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(PIQ_out) <- c("TP","FP","FN","TN","i")

    for(i in Peak_sequence){
        sub <- peak_regions[which(peak_regions$score >=i)]

        sub <- intersect_with_metadata(sub,PIQ_bed)

        sub$Motif <- motif_match_mat[names(sub),TF]
        sub2 <- sub[which(sub$Motif==1)]
        
        int <- intersect_with_metadata(sub2,TF_bed)
        #print(int)
        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub2),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub2),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        PIQ_out <- rbind(PIQ_out,cbind(TP,FP,FN,TN,i))
        print(paste0("Finished ", grep(i, Peak_sequence), " out of ",length(Peak_sequence)))
    }
    colnames(PIQ_out) <- c("TP","FP","FN","TN","Cutoff")
    index <- which(PIQ_out$TP == 0)
    PIQ_out$Precision <- PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FP)
    PIQ_out$Precision[index] <- 0
    PIQ_out$Recall <-  PIQ_out$TP/( PIQ_out$TP+ PIQ_out$FN)
        
    PIQ_out$Baseline <- length(TF_bed)/length(peak_regions)
    PIQ_out$Area <- abs(trapz(PIQ_out$Recall, PIQ_out$Precision))
    PIQ_out$TF <- TF
    PIQ_out$Method <- "HINT"

  ##  str(PIQ_out)

    print("Finished HINT")

    motif_match_out <- data.frame(stringsAsFactors=F)
    motif_match_out <- rbind(motif_match_out,cbind(length(TF_bed),Inf,0,(length(peak_regions)-length(TF_bed)),0))
    colnames(motif_match_out) <- c("TP","FP","FN","TN","i")
    
    for(i in Peak_sequence){
        sub <- peak_regions[which(peak_regions$score >=i)]
        sub$Motif <- motif_match_mat[names(sub),TF]
        sub2 <- sub[which(sub$Motif==1)]
        
        int <- unique(intersect_with_metadata(sub2,TF_bed))

        TP <- length(unique(names(int)))
        FP <- length(unique(setdiff(names(sub2),names(TF_bed))))
        FN <- length(unique(setdiff(names(sub2),x=names(TF_bed))))
        TN <- length(peak_regions)-(TP+FP+FN)
        motif_match_out <- rbind(motif_match_out, cbind(TP,FP,FN,TN,i))
        print(paste0("Finished ", grep(i, Peak_sequence), " out of ",length(Peak_sequence)))
    }
    colnames(motif_match_out) <- c("TP","FP","FN","TN","Cutoff")

    index <- which(motif_match_out$TP == 0)
    motif_match_out$Precision <- motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FP)
    motif_match_out$Precision[index] <- 0
    motif_match_out$Recall <-  motif_match_out$TP/( motif_match_out$TP+ motif_match_out$FN)
    
        
    motif_match_out$Baseline <- length(TF_bed)/length(peak_regions)
    motif_match_out$Area <- abs(trapz(motif_match_out$Recall, motif_match_out$Precision))
    motif_match_out$TF <- TF
    motif_match_out$Method <- "Motifmatchr"

##    str(motif_match_out)
    print("Finished motif match")
    print(paste0("MotifMatch Area for ",TF,":",round(unique(motif_match_out$Area),2)))
    print(paste0("HINT Area for ",TF,":",round(unique(PIQ_out$Area),2)))

    PIQ_out$index <- grep(-0.5,PIQ_out$Cutoff)
    motif_match_out$index <- grep(-0.5, motif_match_out$Cutoff)+nrow(PIQ_out)
     
    final_out <- rbind(PIQ_out,motif_match_out)
    final_out$Difference <- final_out$Area-final_out$Baseline
return(final_out)
}



jaspar_pfm_to_pwm <- function(pfmin,base=c(2,10,"natural"),pseudocounts=0.01,base_prob=0.25) {
    pwmnorm=t(t(pfmin)/colSums(pfmin))
    if(base==2){
        pwmnorm = log2(pwmnorm+pseudocounts)-log2(base_prob)
    } else if(base==10){
        pwmnorm = log10(pwmnorm+pseudocounts)-log10(base_prob)
    } else if(base=="natural"){
        pwmnorm=log(pwmnorm+pseudocounts)-log(base_prob)}

   
return(pwmnorm)
}

CISBP_PFM_to_PWM <- function(PFM,PFM_ID,metadata_file=NULL,outfile=NULL,base=2,base_prob=0.25, pseudocounts=0.01) {
    require(stringr)
    if(is.null(metadata_file)){
        PWM_name <- str_extract(PFM,"M[0-9][0-9][0-9][0-9]_")
        PWM_name <- gsub("\\.txt","",PWM_name)
    } else{ metadata_df <- read.delim(metadata_file,header=T, stringsAsFactors=F)
##            str(metadata_df)
            ID <- str_extract(PFM,"M[0-9][0-9][0-9][0-9]_")
            ID <- paste0(ID, unlist(strsplit(PFM,ID))[2])
            ##print(ID)
            ID <- gsub("\\.txt","",ID)
            ##print(ID)
            PWM_name <- metadata_df[grep(ID, metadata_df$Motif_ID),"Motif_Name"]
##            print(PWM_name)
           }                     

    PFM <- t(read.table(PFM, header=T)[,-1])
    PWM_mat <- jaspar_pfm_to_pwm(PFM,base,pseudocounts=pseudocounts,base_prob=base_prob)
    
    if(is.null(outfile)){
        PWM_list <- list()
        if(length(PWM_name) >0){
            for(i in PWM_name){
                PWM_list[[i]] <- PWM_mat}
            return(PWM_list)}
        else { PWM_list <- list(PWM_mat)
                                                    return(PWM_list)}
    } else{
        PWM_out <- list(paste0(">",PWM_name))
        for(i in 1:nrow(PWM_mat)){ vec <- paste0(rownames(PWM_mat)[i]," \t"," [ ", paste0(unlist(PWM_mat[i,]),collapse=" ")," ] ")
                                   PWM_out <- c(PWM_out,list(vec))}
        writeLines(PWM_out, outfile)
        print("Writing PWM to file")
    }}

Jaspar_PWM_to_file <- function(PWM, PWM_name, outfile){
    PWM <- as.matrix(PWM)
        PWM_out <- list(paste0(">",PWM_name))
        for(i in 1:nrow(PWM)){ vec <- paste0(rownames(PWM)[i]," \t"," [ ", paste0(unlist(PWM[i,]),collapse=" ")," ] ")
                               PWM_out <- c(PWM_out,list(vec))}
    if(is.null(outfile) ==FALSE){
        writeLines(PWM_out, outfile)
        print("Writing PWM to file")} else{ return(PWM_out)}}


get_direct_druggable_candidates <- function(Achilles_df, NTP_df,DGIDB_df, method=c("zscore","cutoff"), cutoff=NULL,n_candidates=100){
    require(dplyr)
    print("grouping")

    NTP_df <- NTP_df[which(NTP_df[,1] %in% colnames(Achilles_df)),]
  ##  print(nrow(Achilles_df))
    Achilles_df <- Achilles_df[intersect(rownames(Achilles_df),DGIDB_df[,1]),]

  ##  print(nrow(Achilles_df))

    grouped_df <- as.data.frame(do.call("cbind", lapply(unique(NTP_df$Prediction), function(x) rowMedians(Achilles_df[NTP_df[which(NTP_df$Prediction==x),1]]))))
    rownames(grouped_df) <- rownames(Achilles_df)
    colnames(grouped_df) <- unique(NTP_df$Prediction)

    if(method=="cutoff"){
        cutoff <- ifelse(is.null(cutoff),-0.5,cutoff)
        print(paste0("Using fixed essentiality cutoff of ",cutoff))
        
        grouped_df_list <- lapply(colnames(grouped_df), function(x) intersect(rownames(sort_df(grouped_df,x,sort_order=FALSE)), rownames(grouped_df[which(grouped_df[,x] <=cutoff),]))[1:n_candidates])
        names(grouped_df_list) <- colnames(grouped_df)
        str(grouped_df_list)
        final_df <- do.call("rbind",lapply(names(grouped_df_list), function(x) data.frame(grouped_df_list[[x]],grouped_df[grouped_df_list[[x]],x],x,stringsAsFactors=F)))
        colnames(final_df) <- c("Gene","Essentiality_in_cluster","Cancer_Type")
           
    } else if(method=="zscore"){
        grouped_df_zscore <- as.data.frame(zscore_df(grouped_df,"row"))
        grouped_df_list <- lapply(colnames(grouped_df), function(x) intersect(rownames(sort_df(grouped_df,x,FALSE))[1:5000],x=rownames(sort_df(grouped_df_zscore,x,FALSE))[1:1000])[1:n_candidates])
        names(grouped_df_list) <- colnames(grouped_df)
##        str(grouped_df_list)
        final_df <- do.call("rbind",lapply(names(grouped_df_list), function(x) data.frame(grouped_df_list[[x]],grouped_df[grouped_df_list[[x]],x],x,stringsAsFactors=F)))
        colnames(final_df) <- c("Gene","Essentiality_in_cluster","Cancer_Type")
        
    }
    final_df$Method <- method
    final_df$Candidates <- paste0(final_df$Cancer_Type,"~",final_df$Gene)
    final_df <-  final_df %>% group_by(Cancer_Type) %>% mutate(Final_Rank=rank(Essentiality_in_cluster)) %>% data.frame
    colnames(final_df)[1] <- "Gene2"

     return(final_df)
}


nawaf_motif_match <- function(PFM_list, peaks,genome=c("hg38","hg19"),output="positions"){
    require(motifmatchr)
    out <- unlist(matchMotifs(PFM_list,peaks,genome=genome,out=output))
    out$Motif <- names(out);
    out$TF <- unlist(lapply(out$Motif, function(x) unlist(strsplit(x,"_"))[1]))
    return(out)}


plot_NTP_heatmaps <- function(NTP_output,metadata_df,column="Disease",pals="Greys",dir=1,title=NULL,is_patient=FALSE,drop_empty=FALSE){

    require(dplyr)
    require(ggplot2)
    df <- as.data.frame(table(NTP_output[c("Prediction","Subtype","Signif")]),stringsAsFactors=F)
    colnames(df)[2] <- "Cancer"
    df <- df %>% group_by(Cancer) %>% mutate(Freq_Norm=Freq/sum(Freq)) %>% data.frame
    missing_preds <- setdiff(metadata_df[,column],unique(df$Prediction))

    if(is_patient==TRUE){
        missing_preds_df <- as.data.frame(do.call("rbind", lapply(missing_preds, function(x) data.frame(x, unique(metadata_df[,column]),TRUE,0,0))))} else {   missing_preds_df <- as.data.frame(do.call("rbind", lapply(missing_preds, function(x) data.frame(x, unique(metadata_df$Subtype),TRUE,0,0))))}
    
    if(nrow(missing_preds_df) >0){
        colnames(missing_preds_df) <- c("Prediction","Cancer","Signif","Freq","Freq_Norm")
        missing_preds_df[1:2] <- apply(missing_preds_df[1:2],2,as.character)
        df <- rbind(df,missing_preds_df,stringsAsFactors=F)} else { df <- df}
    df[4:5] <- apply(df[4:5],2,as.numeric)
    str(df)

    


   centroid_df <- as.data.frame(do.call("rbind",lapply(unique(metadata_df[,column]), function(x) colMeans(dplyr::filter(metadata_df, get(column)==x)[,1:2]))),stringsAsFactors=F)
    centroid_df[,column] <- unique(metadata_df[,column])

    centroid_df <- centroid_df[hclust(dist(centroid_df[,1:2]))$order,]
    factor_levels <- centroid_df[,column]
    df$Prediction <- factor(df$Prediction, levels=factor_levels)
    

###Quick hack of factor levels for better plots #####

    df_mat <- matrix(0, nrow=length(unique(df$Prediction)), ncol=length(unique(df$Cancer)))
    rownames(df_mat) <- unique(df$Prediction)
    colnames(df_mat) <- unique(df$Cancer)
    index <- which(df[,5] !=0) 
    for(n in index){
        i=as.character(df[n,1])
        j=df[n,2]
        val=df[n,5]
        df_mat[i,j] <- val
    }


    df_mat <- df_mat[as.character(factor_levels),]

    j_blacklist <- c()
    indices <- c()
    for(i in 1:nrow(df_mat)){
        best <- sort(df_mat[i,],sort_order=TRUE)
        best_available <- best[setdiff(names(best),j_blacklist)]
        if(best_available[1] >=0.2) {j_blacklist <- c(j_blacklist,names(best_available[1])); indices <- c(indices,rownames(df_mat)[i])} else{ j_blacklist <- j_blacklist; indices <- indices}
        
            

  ##      print(j_blacklist)
    }

##    print(indices)
    for(i in setdiff(unique(df$Cancer),j_blacklist)){
##        print(i)
        best <- sort(df_mat[,i],sort_order=TRUE)
##        #print(best)
        if(best[1] !=0){ insert_loc <- grep(names(best)[1], indices)
                         str(insert_loc)
                         j_blacklist <- c(j_blacklist[1:insert_loc], i, j_blacklist[(insert_loc+1):length(j_blacklist)])} else{ j_blacklist <- c(j_blacklist,i)
                                                                                                                            }
##        print(j_blacklist)
    }

    
##    factor_levels2 <- c("LIHC", "TGCT", "STAD", "CHOL", "COAD", "BRCA", "SKCM", "THCA", "BLCA", "UCEC", "LUAD", "ESCA", "LUSC", "HNSC", "CESC", "PRAD", "KIRC", "KIRP", "ACC","MESO", "PCPG","LGG")
    factor_levels2 <- j_blacklist
   ## print(df_mat[factor_levels[1:10],factor_levels2[1:10]])
   ## print(df_mat[factor_levels[1:10],j_blacklist[1:10]])

    df$Cancer <- factor(as.character(df$Cancer), levels=factor_levels2)
    print("Finished arranging heatmaps")

    df2 <- dplyr::filter(df, Freq!=0)

    if(drop_empty==TRUE){
        df <- dplyr::filter(df, Freq !=0)
    print("Dropping missing values")} else{ df <- df}

    str(df)

    p <- ggplot(df, aes(y=Prediction,Cancer,fill=Freq_Norm))+geom_tile()+scale_fill_distiller(palette=pals,direction=dir)+geom_point(data=df2,aes(size=Freq),fill="white",color="red",shape=1)+theme_classic()+theme(panel.background=element_rect(fill="gray96"),axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))
    if(is.null(title)){ p <- p} else{ p <- p+labs(title=title)}
    return(p)
    }


    
clean_DGIDB <- function(DGIDB) {
    numeric_index <- which(is.na(as.numeric(DGIDB$drug_name))==FALSE)
    chembl_index <- grep("CHEMBL", DGIDB$drug_name)


    DGIDB[chembl_index,"drug_name"] <- DGIDB[chembl_index,"drug_claim_primary_name"]
    DGIDB[numeric_index,"drug_name"] <- DGIDB[numeric_index,"drug_claim_primary_name"]
 

    return(DGIDB)
}


plot_motif_match_comparisons <- function(validation_output,cutoff=-0.5){
    require(dplyr)
    require(ggplot2)
    require(ggsignif)

##    if(sort(unique(validation_output$Method))==c("HINT","Motifmatchr")){ meta_col <- "Footprint_Depth"}


    df <- validation_output
    df <- set_highlight(df,"Cutoff",cutoff,"Size",2.5,0.1)
    sub_df <- unique(df[8:ncol(df)])
##    str(sub_df)
    sub_df <- unique(sub_df[setdiff(colnames(sub_df),c("Size","index"))])
 ##   str(sub_df)

    comparisons <- combn(unique(validation_output$Method),2,simplify=FALSE)

    pr_df <- dplyr::filter(df, Size==2.5)
    ##   str(pr_df)

    pr_df$Precision_Ratio <- pr_df$Precision/pr_df$Baseline

    
    

    p <- ggplot(df, aes(Recall,Precision,fill=Method, alpha=0.4,color=Method))+geom_line()+geom_point(size=df$Size)+facet_wrap(TF~.)+theme(axis.text.x=element_text(face="bold",size=6),axis.text.y=element_text(face="bold",size=6),strip.text = element_text(colour = "black", face = "bold",size=10))+geom_hline(aes(yintercept=Baseline,alpha=0.4),color="gray50",linetype=2)+facet_wrap(TF~.)
    p2 <- ggplot(sub_df, aes(Method, Area,label=TF,fill=Method))+geom_point(aes(color=Method,alpha=0.4))+geom_hline(aes(yintercept=0),color="red", linetype=2, size=0.5,alpha=0.4)+geom_boxplot(width=0.05,fill="white",alpha=0.1)+geom_violin(aes(alpha=0.2),trim=FALSE)+geom_text_repel(size=2,fontface="bold")+geom_signif(test="wilcox.test",comparisons=comparisons,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25),step_increase=0.1,tip_length=0,margin_top=0.4)+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
    p3 <- ggplot(sub_df, aes(Method, Difference,label=TF,fill=Method))+geom_point(aes(color=Method,alpha=0.4))+geom_hline(aes(yintercept=0),color="red", linetype=2, size=0.5,alpha=0.4)+geom_boxplot(width=0.05,fill="white",alpha=0.1)+geom_violin(aes(alpha=0.2),trim=FALSE)+geom_text_repel(size=2,fontface="bold")+geom_signif(test="wilcox.test",comparisons=comparisons,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25),step_increase=0.1,tip_length=0,margin_top=0.4)+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))

    p4 <- ggplot(pr_df, aes(Method, Precision,label=TF,fill=Method))+geom_point(aes(color=Method,alpha=0.4))+geom_hline(aes(yintercept=0),color="red", linetype=2, size=0.5,alpha=0.4)+geom_boxplot(width=0.05,fill="white",alpha=0.1)+geom_violin(aes(alpha=0.2),trim=FALSE)+geom_text_repel(size=2,fontface="bold")+geom_signif(test="wilcox.test",comparisons=comparisons,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25),step_increase=0.1,tip_length=0,margin_top=0.4)+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
    p5 <- ggplot(pr_df, aes(Method, Precision_Ratio,label=TF,fill=Method))+geom_point(aes(color=Method,alpha=0.4))+geom_hline(aes(yintercept=0),color="red", linetype=2, size=0.5,alpha=0.4)+geom_boxplot(width=0.05,fill="white",alpha=0.1)+geom_violin(aes(alpha=0.2),trim=FALSE)+geom_text_repel(size=2,fontface="bold")+geom_signif(test="wilcox.test",comparisons=comparisons,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25),step_increase=0.1,tip_length=0,margin_top=0.4)+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
    p6 <- ggplot(pr_df, aes(Method, Recall,label=TF,fill=Method))+geom_point(aes(color=Method,alpha=0.4))+geom_hline(aes(yintercept=0),color="red", linetype=2, size=0.5,alpha=0.4)+geom_boxplot(width=0.05,fill="white",alpha=0.1)+geom_violin(aes(alpha=0.2),trim=FALSE)+geom_text_repel(size=2,fontface="bold")+geom_signif(test="wilcox.test",comparisons=comparisons,map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05,"."=0.1,"_"=0.25),step_increase=0.1,tip_length=0,margin_top=0.4)+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
    

    print(p)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    print(p6)
}


jaccard <- function(vec1,vec2){
    out <- length(intersect(vec1,vec2))/length(union(vec1,vec2))
    return(out)}

gr_jaccard <- function(gr1, gr2,method=c("standard","bedtools")){

    if(method =="bedtools"){
    int <- sum(width(reduce(intersect(gr1,gr2))))
    combined <- sum(width(reduce(union(gr1,gr2))))
    out <- int/(combined-int)} else if (method=="standard"){
        int <- sum(width(reduce(intersect_with_metadata(gr1,gr2))))
        combined <- sum(width(reduce(union(gr1,gr2))))
        out <- int/combined}
    

##    print(int)
##    print(combined)

    return(out)}

    
    
plot_NTP_comparison <- function(NTP1,NTP2,metadata_df,label1="Method1",label2="Method2"){
    require(ggplot2)
    require(reshape2)

    

    NTP1 <- dplyr::filter(NTP1, Signif==TRUE)
    NTP2 <- dplyr::filter(NTP2, Signif==TRUE)


    NTP1_back <- NTP1
    NTP2_back <- NTP2

    assignment_df <- merge(as.data.frame(table(NTP1$Prediction)),as.data.frame(table(NTP2$Prediction)),by="Var1",all=TRUE)
    assignment_df[,1] <- as.character(assignment_df[,1])
    assignment_df$Diff <- assignment_df[,3]-assignment_df[,2]
    str(assignment_df)
    rownames(assignment_df) <- assignment_df[,1]

    NTP1 <- NTP1[c("Cell_Line","Prediction","Subtype","Detailed_Cancer")]
    NTP2 <- NTP2[c("Cell_Line","Prediction","Subtype","Detailed_Cancer")]

    NTP_comparison <- merge(NTP1, NTP2, by="Cell_Line",all=TRUE)

    NTP_comparison[,2][is.na(NTP_comparison[,2])] <- "Unclassified"
    NTP_comparison[,5][is.na(NTP_comparison[,5])] <- "Unclassified"
    str(NTP_comparison)
    NTP_comparison_table <- as.data.frame(table(NTP_comparison[c("Prediction.x","Prediction.y")]))
    colnames(NTP_comparison_table)[1:2] <- c(label1,label2)
    NTP_comparison_table[1:2] <- apply(NTP_comparison_table[1:2],2,as.character)
    NTP_comparison_table <- NTP_comparison_table %>% group_by(get(label2)) %>% mutate(Freq_Norm=Freq/sum(Freq)) %>% data.frame
    
    ##NTP_comparison_table <- NTP_comparison_table %>% group_by(get(label2)) %>% mutate(Diff=Freq-nrow(dplyr::filter(NTP1, Prediction==get(label1)))) %>% data.frame
    str(NTP_comparison_table)

    index <- which(NTP_comparison_table$Old == NTP_comparison_table$New)

    NTP_comparison_table$Diff <- NA
    NTP_comparison_table[index,"Diff"] <- assignment_df[NTP_comparison_table[index,"Old"],"Diff"]

    print(unique(NTP_comparison_table[,1]))
    print(unique(NTP_comparison_table[,2]))
    df2 <- dplyr::filter(NTP_comparison_table, !is.na(Diff))
    print(df2)
    print(table(NTP_comparison_table$Diff))


    p <- ggplot(NTP_comparison_table, aes(get(label1),get(label2),fill=Freq_Norm))+geom_tile()+scale_fill_distiller(palette="Greys",direction=dir)+theme_classic()+theme(panel.background=element_rect(fill="gray96"),axis.text.x=element_text(face='bold',size=10,angle=30,hjust=0.95),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=16))+geom_point(data=df2,aes(size=Diff),fill="white",color="red",shape=1)+xlab(label1)+ylab(label2)
    

    print(p)
    p2 <- plot_NTP_heatmaps(NTP1_back,metadata_df)+labs(title=label1)
    p3 <- plot_NTP_heatmaps(NTP2_back,metadata_df)+labs(title=label2)


    print(p2)
    print(p3)
}

plot_carnets_candidate_comparison <- function(carnets1, carnets2,known_candidate_list,label1,label2,max_rank=250,plot_mode=c("all","none","some"),fontsize=10,rank_col1="Final_Rank",rank_col2="Final_Rank"){
    require(dplyr)
    require(reshape)
    require(ggplot2)

    carnets1$Candidate_status <- carnets1$Candidates %in% known_candidate_list
    carnets2$Candidate_status <- carnets2$Candidates %in% known_candidate_list
    print("merging tables")

    required_cols1 <- c("Gene2","Cancer_Type","Candidates",rank_col1,"Candidate_status")
    required_cols2 <- c("Gene2","Cancer_Type","Candidates",rank_col2,"Candidate_status")


    carnets1 <- unique(carnets1[required_cols1])
    carnets2 <- unique(carnets2[required_cols2])
    colnames(carnets1)[4] <- colnames(carnets2)[4] <- "Final_Rank"
    carnets1 <- dplyr::filter(carnets1,!Cancer_Type=="",Final_Rank <max_rank)

    carnets2 <- dplyr::filter(carnets2,!Cancer_Type=="", Final_Rank <max_rank)

    carnets_final <- merge(carnets1,carnets2, by="Candidates",all=TRUE)
    print("plotting")

    if(plot_mode=="none"){ carnets_final <- carnets_final} else if( plot_mode=="all"){

    carnets_final$Final_Rank.x[is.na(carnets_final["Final_Rank.x"])] <- max_rank
    carnets_final$Final_Rank.y[is.na(carnets_final["Final_Rank.y"])] <- max_rank} else if(plot_mode=="some"){
        index_x <- intersect(which(carnets_final["Final_Rank.x"] <=10),which(is.na(carnets_final["Final_Rank.y"])==TRUE))
        index_y <- intersect(which(carnets_final["Final_Rank.y"] <=10),which(is.na(carnets_final["Final_Rank.x"])==TRUE))
        carnets_final[index_x,"Final_Rank.y"] <- max_rank
        carnets_final[index_y,"Final_Rank.x"] <- max_rank
    }

    carnets_final$Candidate_status <- carnets_final$Candidates %in% known_candidate_list
    carnets_final$Cancer_Type <- unlist(lapply(1:nrow(carnets_final), function(x) unique(as.character(na.omit(c(carnets_final[x,"Cancer_Type.x"],carnets_final[x,"Cancer_Type.y"]))))))
    carnets_final$Gene2 <- unlist(lapply(carnets_final$Candidates, function(x) unlist(strsplit(x,"~"))[3]))
##    str(carnets_final)
    sub_carnets <- carnets_final[unique(c(which(carnets_final["Final_Rank.x"] <=10),which(carnets_final["Final_Rank.y"] <=10))),]
    sub_carnets$col <- "black"
    sub_carnets$fontface <- "plain"
    sub_carnets[which(sub_carnets$Candidate_status==TRUE),"col"] <- "red"
    sub_carnets[which(sub_carnets$Candidate_status==TRUE),"fontface"] <- "bold"
    ##    str(sub_carnets)



    p <- ggplot(carnets_final, aes(Final_Rank.x,Final_Rank.y, label=Gene2))+geom_point(size=2,alpha=0.5,color="gray50",aes(shape=Candidate_status,fill=Cancer_Type))+facet_wrap(Cancer_Type~.,ncol=4,scales="free")+geom_text_repel(data=sub_carnets,size=3, fontface=sub_carnets$fontface,show.legend=F,color=sub_carnets$col,segment.size=0.1,segment.alpha=0.3)+scale_TF_score("~/ANF_Final_Cluster_key.txt",extension=".txt")+theme_bw()+theme(legend.position="right",axis.text.x=element_text(face='bold',size=fontsize),axis.text.y=element_text(face='bold',size=fontsize),legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=8,face="bold"),strip.text.x = element_text(colour = "black", face = "bold",size=8),strip.text.y = element_text(colour = "black", face = "bold",size=8))+guides(fill=guide_legend(ncol=1),color=guide_legend(ncol=1))+scale_x_reverse()+scale_y_reverse()+xlab(label1)+ylab(label2)+scale_shape_manual(values=c(21,22))

    print(p)
}

get_potential_cell_line_controls <- function(carnets_df,NTP_df,Achilles_df,cancer_subset=NULL,max_rank=3){

    check <- grep("GEO", colnames(NTP_df),value=TRUE)
##    str(check)
    if(length(check) ==0){ print("No GEO data in the NTP file. Can't prioritize cell lines")} else{
       ## NTP_df <- dplyr::filter(NTP_df, nchar(get(check)) !=0)
##        print(dplyr::filter(NTP_df, grepl("LUAD-PCPG",Prediction)))
        NTP_df <- do.call("rbind", lapply(unique(NTP_df$Cancer_Type), function(x) sort_df(dplyr::filter(NTP_df, Cancer_Type==x),check,TRUE)))
        NTP_df_back <- NTP_df

        NTP_df[,check] <- as.numeric(NTP_df[,check])
        NTP_df <- dplyr::filter(NTP_df, Cell_Line %in% colnames(Achilles_df))
        NTP_df <- NTP_df %>% group_by(Cancer_Type) %>% mutate(Cell_Line_Rank=(length(get(check))+1)-rank(get(check))) %>% data.frame
        
    }

    
    
    if(is.null(cancer_subset) == TRUE){
        carnets_df <- dplyr::filter(carnets_df, Final_Rank <=max_rank)
        NTP_df <- dplyr::filter(NTP_df, Signif==TRUE,Cell_Line_Rank <=max_rank)
    } else{ carnets_df <- dplyr::filter(carnets_df, Final_Rank <=max_rank,Cancer_Type %in% cancer_subset)
            NTP_df <- dplyr::filter(NTP_df, Signif==TRUE,Cancer_Type %in% cancer_subset, Cell_Line_Rank <=max_rank)
            
        }

    print(NTP_df)

    rownames(NTP_df) <- NTP_df[,1]

##    str(NTP_df)
##    str(carnets_df)

    controls_df <- expand.grid(NTP_df[,1], carnets_df[,2],2, stringsAsFactors=F)[1:2]
    controls_df$Cell_Line_Cluster <- NTP_df[controls_df[,1],"Cancer_Type"]
    controls_df$Gene_Cluster <- unlist(lapply(controls_df[,2], function(x) paste0(unique(dplyr::filter(carnets_df, Gene2==x)$Cancer_Type),collapse="_")))

    print(table(controls_df$Gene_Cluster))
    print(table(controls_df[,2] %in% rownames(Achilles_df)))

    controls_df$Essentiality <- NA
    controls_df$Essentiality_in_target_cluster <- NA

##    str(controls_df)
##    str(NTP_df)
    colnames(controls_df)[1:2] <- c("Cell_Line","Gene")

##    print(controls_df[577:580,])

    for(i in 1:nrow(controls_df)){ if(controls_df[i,1] %in% colnames(Achilles_df)){

                                      ##print(paste0("Working ",i))
                                      in_cluster <- dplyr::filter(NTP_df_back, Cancer_Type==controls_df[i,4], Cell_Line %in% colnames(Achilles_df))[,1]
                                     ## print(in_cluster)
                                      controls_df[i,"Essentiality"] <- Achilles_df[controls_df[i,2],controls_df[i,1]]
                                      controls_df[i,"Essentiality_in_target_cluster"] <- median(unlist(Achilles_df[controls_df[i,2],in_cluster]),na.rm=TRUE)


                                  } else {

                                          controls_df[i,"Essentiality"] <- median(unlist(Achilles_df[controls_df[i,2],dplyr::filter(NTP_df_back, Cancer_Type==controls_df[i,3],Cell_Line %in% colnames(Achilles_df))[,1]]),na.rm=TRUE)
                                          in_cluster <- dplyr::filter(NTP_df_back, Cancer_Type==controls_df[i,4], Cell_Line %in% colnames(Achilles_df))[,1]
                                        ##  print(in_cluster)
                                          
                                          controls_df[i,"Essentiality_in_target_cluster"] <- median(unlist(Achilles_df[controls_df[i,2],in_cluster]),na.rm=TRUE)
                                          
                                  }
#print(i)
                               }

    controls_df$Differential_Essentiality <- controls_df$Essentiality-controls_df$Essentiality_in_target_cluster
    controls_df$Essentiality_Ratio <- controls_df$Essentiality_in_target_cluster/controls_df$Essentiality
    controls_df$Good_Control <- FALSE
    controls_df[intersect(which(controls_df$Essentiality_Ratio >= 2), which(controls_df$Differential_Essentiality >=0.1)),"Good_Control"] <- TRUE

##    controls_df$Essentiality <- Achilles_df[controls_df[,2],controls_df[,1]]
##    controls_df$Essenotiality_target_cluster <- Achilles_df[controls_df[,2],unlist(lapply(controls_df[,1], function(x) dplyr::filter(NTP_df,Prediction==x)))
    return(controls_df)
    
}

get_promoters <- function(gene_gr,upstream=2500,downstream=100){
    require(dplyr)
    require(reshape2)
    source("/home/forbesa1/Andre_F_functions.R")

    coords <- data.frame(as.character(gene_gr@seqnames),gene_gr@ranges@start, gene_gr@ranges@start+gene_gr@ranges@width,as.character(gene_gr@strand),stringsAsFactors=F)
    colnames(coords) <- c("chr","start","end","strand")
    coords <- cbind(coords, metadata <- as.data.frame(gene_gr@elementMetadata@listData,stringsAsFactors=F))
##    str(coords)
    coords1 <- dplyr::filter(coords,strand=="+")
    coords2 <- dplyr::filter(coords,strand=="-")


    coords3 <- dplyr::filter(coords,!strand %in% c("+","-"))


    coords1$end <- coords1$start+downstream
    coords1$start <- coords1$start-upstream
    
    coords2$start <- coords2$end-downstream
    coords2$end <- coords2$end+upstream

    coords3$end <- coords3$start+downstream
    coords3$start <- coords3$start-upstream

    out <- rbind(coords1,coords2,coords3,stringsAsFactors=F)
    out <- table_to_granges(out,strand=4)

    return(out)


}
                         

HiCPro_to_granges <- function(HiCPro_file,filter_trans=TRUE,buffer=50,ranges_1=1:4,ranges_2=5:8){
    gz_check <- grep("\\.gz", HiCPro_file)
    if(length(gz_check) ==0){ print(paste0("Working with ",HiCPro_file))

                              system(paste0("awk -v OFS='\t\' '{print $1, $2-",buffer,",$2+",buffer,",$3",",$4,$5-",buffer,",$5+",buffer,",$6,$7,$8}' ",HiCPro_file," > HiCPro_bedpe.bedpe"))


                          } else {print("Unzipping file, hope there's enough storage")
                                  system(paste0("zcat ", HiCPro_file," | awk '{OFS=\"\t\"} {print $2,$3,$4,$5,$6,$7,$8,$1}' > intermediate_HiCPro_file.txt"))
                                  system(paste0("awk -v OFS='\t' '{print $1, $2-",buffer,",$2+",buffer,",$3","$4,$5-",buffer,",$5+",buffer,",$6,$7,$8}' intermediate_HiCPro_file.txt > HiCPro_bedpe.bedpe"))

                                  system("head HiCPro_bedpe.bedpe")
                              }

    links <- bedpe_to_granges("HiCPro_bedpe.bedpe", ID_col=10,ranges_1,ranges_2)
    print(links)


}

parse_HiCPro <- function(HiCPro_bedpe, TCGA_peaks, promoter_gr){

    require(dplyr)

    sub <- unlist(HiCPro_bedpe)
    sub$Target <- "Unknown"
    sub$Enhancer <- "Unknown"
    sub$Element <- "Unknown"

    enhancer_gr <- setdiff_with_metadata(TCGA_peaks, promoter_gr)
##    print(enhancer_gr)


    sub_e <- intersect_with_metadata(sub,enhancer_gr)
    enhancer_hits <- as.data.frame(findOverlaps(sub_e,enhancer_gr))
    sub_p <- intersect_with_metadata(sub,promoter_gr)
    promoter_hits <- as.data.frame(findOverlaps(sub_p,promoter_gr))

    ##str(enhancer_hits)
    ##str(promoter_hits)

    sub_e[enhancer_hits[,1]]$Enhancer <- enhancer_gr[enhancer_hits[,2]]$name
    sub_e$Element <- "Enhancer"
    sub_p[promoter_hits[,1]]$Target <- promoter_gr[promoter_hits[,2]]$HUGO
    sub_p$Element <- "Promoter"

    int <- c(sub_p, sub_e)

    int2 <- c(int, setdiff_with_metadata(sub,int))

    ##print(int2[which(int2$Name=="J00118:213:HFVFGBBXX:5:2111:17898:9825")])


    ##out <- split(int2,int2$Name)    
    out <- gr_to_bed(int2, outfile=NULL, metadata=TRUE,header=TRUE)
    out <- out %>% group_by(Name) %>% mutate(Target=max(gsub("Unknown",NA,Target),na.rm=TRUE)) %>% data.frame

##    print(head(out))
##    print(out[which(out$Name=="J00118:213:HFVFGBBXX:5:2111:17898:9825"),])

    out <- table_to_granges(out,"Strand")

    
##    print(out)

    return(out)   
}


parse_HiCPro_v2 <- function(HiCPro_bedpe, TCGA_peaks, promoter_gr,peak_cutoff=-0.5,peak_height_matrix, sample_ID){

    require(dplyr)

    if(sample_ID %in% colnames(peak_height_matrix)){ print(paste0("Working ",sample_ID))

                                                     TCGA_peaks$score <- peak_height_matrix[names(TCGA_peaks),sample_ID]
                                                     TCGA_peaks <- TCGA_peaks[TCGA_peaks$score >=peak_cutoff]} else { print(paste0(sample_ID, " is not in the provided matrix!"))
                                                                                                                      stop()
                                                                                                                  }

    ##    TCGA_peaks <- setdiff_with_metadata(TCGA_peaks,promoter_gr)

    print("Finding peaks in anchors")

    sub <- HiCPro_bedpe
##    print(sub)
    int <- findOverlaps(sub,TCGA_peaks)
    int_pe <- sub[int@from]
   
    peaks <- names(TCGA_peaks[int@to])
    peaks <- unlist(lapply(peaks, function(x) rep(x,2)))

    
    int_pe <- unlist(int_pe)
    int_pe$Peaks <- peaks
    int_pe <- split(int_pe, int_pe$Name)
    str(peaks)

    print("Finding target genes overlapping anchors")
    int <- findOverlaps(int_pe, promoter_gr)
##    str(int)

    out_pe <- int_pe[int@from]

   ## print(out_pe)
   ## print(sort(unique(unlist(out_pe)$Peaks)))
    genes <- promoter_gr[int@to]$HUGO
    str(genes)
    print("Matching genes to anchors..... could take a while")

##    print(length(out_pe))
    gene_vec <- unlist(lapply(1:length(out_pe), function(x) rep(genes[x],length(out_pe[[x]]))))
##    str(gene_vec)
    str(out_pe@unlistData$Peaks)
    out_pe@unlistData$Target <- gene_vec
    out_pe <- unlist(out_pe)
    out_pe$Element <- "Enhancer"
    out_pe$Link <- paste0(out_pe$Peaks,"~",out_pe$Target)
    ##    print(out_pe)
    index <- which(out_pe$Peaks %in% names(intersect_with_metadata(TCGA_peaks, promoter_gr)))
    if(length(index) >=1){
        out_pe[index]$Element <- "Promoter"} else{ print("None of these is a promoter")}
    names(out_pe) <- out_pe$Name
    out_pe <- split(out_pe, out_pe$Name)
    
    

    ##out <- split(int2,int2$Name)    
    
    ##    print(out)

    print("Finished")

    return(out_pe)   
}

quantify_peak_gene_links <- function(P_G_link_df,true_link_file,sample_ID,filter_enhancers=FALSE, enhancer_gr){


    if(is.character(true_link_file)){    truth_set <- read.table(true_link_file, header=F, sep='\t',stringsAsFactors=F)} else { truth_set <- true_link_file}
    if(filter_enhancers ==TRUE){ truth_set <- truth_set[which(truth_set[,3] =="Enhancer"),]} else{ truth_set <- truth_set}
    
    peaks <- unlist(lapply(rownames(P_G_link_df), function(x) unlist(strsplit(x,"~"))[1]))
##    print(dim(P_G_link_df))
    P_G_link_df <- P_G_link_df[which(peaks %in% enhancer_gr$name),]
##    print(dim(P_G_link_df))
    test_set <- rownames(P_G_link_df)[which(P_G_link_df[,sample_ID]>=1)]
    

    ##str(truth_set)
    ##    str(test_set)
    
    TP <- length(intersect(unique(test_set), unique(truth_set[,4])))
    FP <- length(setdiff(unique(test_set), unique(truth_set[,4])))
    FN <- length(setdiff(unique(test_set), x=unique(truth_set[,4])))

    Precision <- TP/(TP+FP)

    out <- data.frame(TP,FP,FN,Precision,sample_ID)
    print(paste0(sample_ID," Precision=",round(Precision,digits=3)))
    print(paste0("Finished ",sample_ID))
    return(out)

}

link_vector_to_links_grl <- function(link_vec, peak_gr, promoter_gr){
    

    require(GenomicRanges)

    peaks <- unlist(lapply(link_vec, function(x) unlist(strsplit(x,"~"))[1]))
    genes <- unlist(lapply(link_vec, function(x) unlist(strsplit(x,"~"))[2]))
    str(peaks)
    str(genes)

    sub_peak_gr <- peak_gr[peaks]
    
    sub_pro_gr <- promoter_gr[unlist(lapply(genes, function(x) which(promoter_gr$HUGO == x)[1]))]
    genes <- sub_pro_gr$HUGO
    str(genes)


    mcols(sub_peak_gr) <- NULL
    mcols(sub_pro_gr) <- NULL

    sub_peak_gr$Name <- peaks
    sub_peak_gr$Element_Type <- "Peak"
    sub_pro_gr$Name <- genes
    sub_pro_gr$Element_Type <- "Gene"

    sub_peak_gr$Link <- link_vec
    sub_pro_gr$Link <- link_vec

    names(sub_peak_gr) <- NULL
    names(sub_pro_gr) <- NULL

    out <- c(sub_peak_gr, sub_pro_gr)
    out <- out[order(out$Link,start(out))]
    strand(out) <- c("+","-")
    out$Name <- paste0(out$Link,"|",out$Element_Type)
    out <- split(out, out$Link)

    return(out)
}

make_P_G_link_matrix <- function(dir,extension="_true_links.txt",gene=NULL){

    file_list <- list.files(dir, extension, full.names=TRUE)

    mat <- data.frame(0)

    if(is.null(gene)){
        for(i in file_list){

            df <- read.table(i, header=F, sep='\t', stringsAsFactors=F,nrows=1e3)

            sample_ID <- unique(df[,5])
            df[,sample_ID] <- 1

            df <- df[c(colnames(df)[4],sample_ID)]
            rownames(df) <- df[,1]
            df[,1] <- NULL
            mat <- vector_smartmerge(mat,df)}} else{
                for(i in file_list){

                    df <- system(paste0("grep -w ", gene," ",i), intern=TRUE)
                    
                    if(length(df) > 0){
                df <- data.frame(do.call("rbind", lapply(df,function(x) unlist(strsplit(x,"\t")))),stringsAsFactors=F)
                
                sample_ID <- unique(df[,5])
                df[,sample_ID] <- 1
                index <- c(colnames(df)[4], sample_ID)
                df <- unique(df[,index])
                ##str(df)
                rownames(df) <- df[,1]
                df[,1] <- NULL
                mat <- vector_smartmerge(mat,df)
                ##str(mat)
                
            } else{
                    df <- data.frame(0)
                    sample_ID <- gsub(extension,"",last(unlist(strsplit(i,"/"))))
                    
                    df[,sample_ID] <- 0
                    df[,1] <- NULL
                    mat <- vector_smartmerge(mat,df)
                }}}
                
                
                                      
                

    mat <- mat[-1,-1]
    mat[is.na(mat)] <- 0
        return(mat)
       

            }

link_comparison_heatmap <- function(true_link_df, test_df,gene,expression_df,metadata_file,metadata_cols,color_scale,filter_links=4,true_method="HiChIP",test_method="DGTAC",promoter_gr=NULL,fontsize=5){
    require(ComplexHeatmap)
    require(circlize)
    true_link_df <- true_link_df[grep(gene, rownames(true_link_df)),intersect(colnames(true_link_df),colnames(test_df))]
    test_df <- test_df[grep(gene, rownames(test_df)),intersect(colnames(true_link_df),colnames(test_df))]

##    str(true_link_df)
##    str(test_df)

    union_links <- union(rownames(true_link_df),rownames(test_df))

    link_metadata <- data.frame(union_links,"Both",stringsAsFactors=F)
    colnames(link_metadata) <- c("Link","Method")
    rownames(link_metadata) <- link_metadata[,1]
    link_metadata[setdiff(rownames(true_link_df),rownames(test_df)),"Method"] <- "True Links Only"
    link_metadata[setdiff(rownames(true_link_df),x=rownames(test_df)),"Method"] <- "Test Links Only"

    if(is.null(promoter_gr)){link_metadata <- link_metadata} else{
        link_metadata$Peak <- unlist(lapply(link_metadata$Link, function(x) unlist(strsplit(x,"~"))[1]))
        link_metadata$Element <- "Enhancer"
        link_metadata[which(link_metadata$Peak %in% promoter_gr$name),"Element"] <- "Promoter"
        ##str(link_metadata)
        link_metadata$Peak <- NULL}
                                                                   


    sample_metadata <- data.frame(colnames(test_df),metadata_file[colnames(test_df),metadata_cols],stringsAsFactors=F)
    
    colnames(sample_metadata) <- c("Sample",metadata_cols)
    rownames(sample_metadata) <- sample_metadata[,1]
    sample_metadata[,gene] <- unlist(expression_df[gene,sample_metadata$Sample])

    ##str(sample_metadata)
    ##print(summary(sample_metadata[,gene]))
    ##str(link_metadata)

    union_df <- matrix(0, nrow=length(union_links), ncol=ncol(test_df))
    rownames(union_df) <- union_links
    colnames(union_df) <- colnames(test_df)

##    str(union_df)


    true_union <- union_df
    vals <- as.matrix(true_link_df[rownames(true_link_df),colnames(true_link_df)])
##    str(vals)
    true_union[rownames(true_link_df),colnames(true_link_df)] <- vals

    test_union <- union_df
    vals2 <- as.matrix(test_df[rownames(test_df),colnames(test_df)])
##    str(vals2)
    test_union[rownames(test_df),colnames(test_df)] <- vals2

    
##    print(true_union[1:10,1:5])
##    print(test_union[1:10,1:5])

    int <- cbind(as.vector(true_union),as.vector(test_union))
    int <- apply(int,1,function(x) paste0(x[1],x[2]))

    ##str(int)


    converter <- list("00"=0,"01"=-1,"10"=1,"11"=2)

    int <- unlist(converter[int])
##    str(int)

    out <- matrix(int, nrow=nrow(union_df), ncol=ncol(union_df))
    rownames(out) <- rownames(union_df); colnames(out) <- colnames(union_df)

    out <- out[which(rowSums(true_union) >=filter_links),]

##    str(link_metadata)
    sample_metadata <- sort_df(sample_metadata,gene,sort_order=TRUE)
    row_anno <- link_metadata[rownames(out),c(2,2:ncol(link_metadata))]
    col_anno <- sample_metadata[,2:ncol(sample_metadata)]
    colnames(col_anno)[grep(gene, colnames(col_anno))] <- "Gene"

    if(is.null(promoter_gr)){    ha <- rowAnnotation(row_anno)} else{
        link_scale=c("black","gray50"); names(link_scale) <- c("Promoter","Enhancer")
        ha <- rowAnnotation(row_anno,col=list(Element=link_scale[unique(row_anno$Element)]))}
    ha2 <- HeatmapAnnotation(col_anno,col=list(Subtype=color_scale[unique(sample_metadata[,2])], Gene=colorRamp2(c(0,max(unlist(sample_metadata[,gene]))),c("white","deeppink"))),which="column")


##    str(sample_metadata)
    out <- out[,sample_metadata$Sample]
##    print(sample_metadata[colnames(out),])

    ht <- Heatmap(out,top_annotation=ha2,col=colorRamp2(c(-1,0,2),c("blue","white","red")),cluster_rows=TRUE, cluster_columns=FALSE,row_names_gp = gpar(font = 2,fontsize=fontsize),column_names_gp = gpar(font = 2,fontsize=4),heatmap_legend_param=list(at=c(-1,0,1,2),labels=c(paste0("Not in ",true_method),"Not in either",paste0("Not in ",test_method),"In Both")))
    ht_out <- ht+ha
    draw(ht_out,column_title=paste0("Peak-Gene link comparison for all links for ",gene), column_title_gp=gpar(fontsize=15))



    out <- matrix(int, nrow=nrow(union_df), ncol=ncol(union_df))
    rownames(out) <- rownames(union_df); colnames(out) <- colnames(union_df)

    out <- out[rownames(test_df),sample_metadata$Sample]
    row_anno <- link_metadata[rownames(out),c(2,2:ncol(link_metadata))]
    if(is.null(promoter_gr)){    ha <- rowAnnotation(row_anno)} else{
        link_scale=c("black","gray50"); names(link_scale) <- c("Promoter","Enhancer")
         ha <- rowAnnotation(row_anno,col=list(Element=link_scale[unique(row_anno$Element)]))}
    
    ht <- Heatmap(out,top_annotation=ha2,col=colorRamp2(c(-1,0,2),c("blue","white","red")),cluster_rows=TRUE, cluster_columns=FALSE,row_names_gp = gpar(font = 2,fontsize=fontsize),column_names_gp = gpar(font = 2,fontsize=4),heatmap_legend_param=list(at=c(-1,0,1,2),labels=c(paste0("Not in ",true_method),"Not in either",paste0("Not in ",test_method),"In Both")))
    ht_out2 <- ht+ha
    draw(ht_out2,column_title=paste0("Peak-Gene link comparison for predicted links for ",gene), column_title_gp=gpar(fontsize=15))



}


plot_link_comparison <- function(link_comparison_df, metadata_df, columns){
    require(ggplot2)
    require(reshape2)
    require(dplyr)

    

    link_comparison_df[,columns] <- metadata_df[as.character(link_comparison_df[,4]),columns]




    link_comparison_df2 <- lapply(unique(link_comparison_df[,"Method"]), function(x) dplyr::filter(link_comparison_df, Method==x))

##    str(link_comparison_df2)

    for(i in 1:length(link_comparison_df2)){ index <- grep("Precision", colnames(link_comparison_df2[[i]])); colnames(link_comparison_df2[[i]])[index] <- paste0(unique(link_comparison_df2[[i]][,"Method"]),"_Precision")}



    link_comparison_df2 <- merge(link_comparison_df2[[1]], link_comparison_df2[[2]], by="sample_ID")
    link_comparison_df2[,columns] <- metadata_df[as.character(link_comparison_df2[,"sample_ID"]),columns]

##    str(link_comparison_df)
   link_comparison_df3 <- link_comparison_df2
    str(link_comparison_df3)
    

    

   p <- ggplot(link_comparison_df, aes(Precision, fill=Method,alpha=0.2))+geom_density()+theme_classic()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
   p2 <-  ggplot(link_comparison_df, aes(TP, fill=Method,alpha=0.2))+geom_density()+theme_classic()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10)); ggplot(link_comparison_df, aes(Disease,Precision, fill=Method,alpha=0.2))+geom_boxplot()+geom_point(aes(fill=Method),color="black",shape=21,position=position_jitterdodge())+theme_classic()+theme(axis.text.x=element_text(face="bold",size=7,angle=30,hjust=0.95),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
   p3 <- ggplot(link_comparison_df, aes(Subtype,Precision, fill=Method,alpha=0.2))+geom_boxplot()+geom_point(aes(fill=Method),color="black",shape=21,position=position_jitterdodge())+theme_classic()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))
   p4 <-  ggplot(link_comparison_df3, aes(Correlation_Precision,DGTAC_Precision, fill=Subtype,alpha=0.4,label=sample_ID))+geom_point(color="black",shape=21,size=5)+theme_classic()+theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10),strip.text = element_text(colour = "black", face = "bold",size=10))+scale_TF_score(color_list="Corces_Hex_Codes_Cancer.txt")+geom_abline(slope=1,intercept=0, color="gray50", alpha=0.2,linetype=5,size=1)

    print(p)
    print(p2)
    print(p3)
    print(p4)

}



get_best_screen_drugs <- function(screen_df, metadata_df, metadata_column,method=c("cutoff","zscore"),cutoff=NULL,zscore_cutoff=-0.5,abs_cutoff=1000,filter_signif=TRUE,sample_margin="row",zscore_method="mean"){

    require(dplyr)
    print("grouping")

    if(sample_margin == "row"){ screen_df <- as.data.frame(t(screen_df))} else{ screen_df <- screen_df}
    if(filter_signif==TRUE){    metadata_df <- dplyr::filter(metadata_df,Signif==TRUE)} else{
        metadata_df <- metadata_df}


    metadata_df <- metadata_df[which(metadata_df[,"Cell_Line"] %in% colnames(screen_df)),]
##    str(metadata_df)
##    str(screen_df)
    screen_df <- unique(screen_df[intersect(colnames(screen_df),metadata_df[,"Cell_Line"])])

##  str(screen_df)  
  ##  print(nrow(Achilles_df))
   

  ##  print(nrow(Achilles_df))

    groups <- sort(unique(metadata_df[,metadata_column]))
##    print(groups)
    grouped_df <- as.data.frame(do.call("cbind", lapply(groups, function(x) rowMeans(screen_df[dplyr::filter(metadata_df, get(metadata_column)==x)[,"Cell_Line"]],na.rm=TRUE))))

                                                                                                                                rownames(grouped_df) <- rownames(screen_df)
    colnames(grouped_df) <- groups


  

    if(method=="cutoff"){
        cutoff <- ifelse(is.null(cutoff),-0.5,cutoff)
        print(paste0("Using fixed essentiality cutoff of ",cutoff))

             
        grouped_df_list <- lapply(colnames(grouped_df), function(x) intersect(rownames(sort_df(grouped_df,x,sort_order=FALSE)), rownames(grouped_df[which(grouped_df[,x] <=cutoff),])))
        names(grouped_df_list) <- colnames(grouped_df)

        
           
    } else if(method=="zscore"){
        grouped_df_zscore <- as.data.frame(zscore_df(grouped_df,"row",zscore_method))

##        str(grouped_df_zscore)

        df_rel <- as.matrix(grouped_df_zscore)
        df_abs <- as.matrix(grouped_df)


        templates_rel <- lapply(colnames(grouped_df), function(x) names(sort(df_rel[which(df_rel[,x]<=zscore_cutoff),x], sort_order=FALSE)))
        templates_abs <- lapply(colnames(grouped_df), function(x) names(sort(df_abs[,x], sort_order=FALSE))[1:abs_cutoff])


        names(templates_abs) <- colnames(grouped_df)
        names(templates_rel) <- colnames(grouped_df)
        grouped_df_list <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))



        names(grouped_df_list) <- colnames(grouped_df)
    }
##        str(grouped_df_list)

     return(grouped_df_list)
}


get_screen_drug_targets <- function(drug_list, compound_metadata){
    require(dplyr)

    print("Getting drug metadata")
    drug_list <- lapply(drug_list, function(x) do.call("rbind", lapply(x, function(y) dplyr::filter(compound_metadata,drug_name == y))))
##    str(drug_list)

    drug_target_list <- list()
    for(i in names(drug_list)){
        print(paste0("Working ",i))
        sub <- drug_list[[i]]
        sub <- dplyr::filter(sub, target!="")


        if(!is.null(nrow(sub))){

            ranks <- 1
            out <- data.frame(stringsAsFactors=F)
            for(j in unique(sub$drug_name)){ int <- dplyr::filter(sub, drug_name==j);vec <- unique(unlist(strsplit(unique(int$target),", "))); vec <- as.data.frame(cbind(vec, ranks),stringsAsFactors=F);vec[,2] <- as.numeric(vec[,2])
                                             ##str(vec)
                                             out <- rbind(out,vec)
                                             ##str(out)
                                             ranks <- ranks+1}

            print(paste0("Finished ",i))
            ##str(out)
            out <- out %>% group_by(vec) %>% mutate(Rank_Min=min(ranks)) %>% data.frame
            out$Rank_Uniq <- 1:nrow(out)                                               
            
            
            drug_target_list[[i]] <- out
        } else{

                                 vec <- "None"
                                 vec <- data.frame(vec,1e5,1e5,stringsAsFactors=F)
                                 drug_target_list[[i]] <- vec}
    }
        
    ##str(drug_target_list)

    drug_target_list <- lapply(drug_target_list, function(x) dplyr::filter(x, vec!="None"))
                                
     
##    drug_target_list <- lapply(drug_list2, function(x) unique(unlist(strsplit(unique(compound_metadata[which(compound_metadata$drug_name %in% unique(unlist(lapply(x, function(y) unlist(strsplit(y,"\\."))[1])))),"target"]),", "))))


##    drug_target_list <- lapply(drug_target_list, function(x) x[!is.na(x)])

    return(drug_target_list)

}


make_stratified_cv_folds <- function(matrix, metadata_df, metadata_col, n_folds){
    metadata_df <- metadata_df[which(metadata_df[,"Sample"] %in% colnames(matrix)),]
    classes <- unique(metadata_df[,metadata_col])

    class_summary <- data.frame(table(metadata_df[,metadata_col]))
    class_summary[,1] <- as.character(class_summary[,1])



    if((n_folds-1) > min(class_summary[,2])){ print(paste0("Dataset is too imbalanced for ", n_folds," folds. Reduce # of folds or increase size of dataset"))}   else{
           training_list <- list()
           test_list <- list()
           for(i in 1:n_folds){
               train_set <- c()
               test_set <- c()
               for(j in classes){
                   sub_meta <- metadata_df[which(metadata_df[metadata_col] == j),]
                   fold_size <- round(nrow(sub_meta)/n_folds)
                   samples <- sample(sub_meta[,"Sample"],fold_size)



                   test_set <- c(test_set,samples)

               }


               training_list[[paste0("Fold_",i)]] <- setdiff(colnames(matrix),test_set)
               test_list[[paste0("Fold_",i)]] <- test_set

               print(paste0("Finished Processing Fold ",i))
           }
    out_list <- list(test_list,training_list)
    names(out_list) <- c("Test","Train")

    return(out_list)}
}

convert_prism_drugs_to_prism_genes <- function(PRISM_df,PRISM_metadata){
    require(reshape2)
    require(dplyr)
    PRISM_metadata <- dplyr::filter(PRISM_metadata,drug_name %in% colnames(PRISM_df))
    
    targets <- sort(unique(PRISM_metadata$target))

    index <- unlist(lapply(targets, nchar))
##    str(index)
    targets <- targets[which(index!=0)]

##    str(targets)

    target_df <- do.call("cbind",lapply(targets, function(x) rowMeans(PRISM_df[dplyr::filter(PRISM_metadata, target ==x)$drug_name],na.rm=TRUE)))
    colnames(target_df) <- targets
    return(target_df)

}

compare_prism_carnets_drugs <- function(carnets_drugs_df, PRISM_df,metadata_df,metadata_column="Prediction",quantile_cutoff=0.75,alternative_opt=c("two.sided","greater","less")){
    require(dplyr)
    print("grouping")

    PRISM_df <- as.data.frame(t(PRISM_df))
    metadata_df <- dplyr::filter(metadata_df,Signif==TRUE)
    metadata_df <- metadata_df[which(metadata_df[,"Cell_Line"] %in% colnames(PRISM_df)),]

    PRISM_df <- unique(PRISM_df[intersect(colnames(PRISM_df),metadata_df[,"Cell_Line"])])

    
  ##  print(nrow(Achilles_df))
   

  ##  print(nrow(Achilles_df))

    groups <- sort(unique(metadata_df[,metadata_column]))
    grouped_df <- as.data.frame(do.call("cbind", lapply(groups, function(x) rowMeans(PRISM_df[dplyr::filter(metadata_df, get(metadata_column)==x)[,"Cell_Line"]],na.rm=TRUE))))

                                                                                                                          rownames(grouped_df) <- rownames(PRISM_df)
    colnames(grouped_df) <- groups
##    print(groups)
    grouped_df <- zscore_df(grouped_df, "row","mean")
    grouped_list <- lapply(colnames(grouped_df), function(x) data.frame(grouped_df[,x],rownames(grouped_df),1-get_quantile(grouped_df[,x]),stringsAsFactors=F,x))
    grouped_list <- lapply(grouped_list, function(x) sort_df(x,1,sort_order=F))
    grouped_list <- do.call("rbind", grouped_list)
    colnames(grouped_list) <- c("Zscore","Drug","Quantile","Cancer_Type")
    grouped_list <- split(grouped_list,grouped_list$Cancer_Type)
    grouped_list <- lapply(grouped_list, function(x) dplyr::filter(x, Quantile>=quantile_cutoff))

##    str(grouped_list)
##    str(carnets_drugs_df)
##    str(rownames(PRISM_df))
    fisher_list <- lapply(names(grouped_list), function(x) make_and_run_fishers_exact_vectors(unique(dplyr::filter(carnets_drugs_df, Cancer_Type==x)$Gene2), rownames(PRISM_df),grouped_list[[x]][,"Drug"],label=x,alternative_opt=alternative_opt))
    names(fisher_list) <- names(grouped_list)
    out <- do.call("rbind", fisher_list)
    rownames(out) <- 1:nrow(out)
    out$Signif <- out$P_value <=0.2
    out$P_value <- round(out$P_value, digits=3)

    
    return(out)
                          
    

}

compare_prism_carnets_targets <- function(carnets_df, PRISM_df,metadata_df,metadata_column="Prediction",compound_metadata_df,quantile_cutoff=0.75,alternative_opt=c("two.sided","greater","less")){
    require(dplyr)
    print("grouping")

    PRISM_df <- as.data.frame(t(PRISM_df))
    metadata_df <- filter(metadata_df,Signif==TRUE)
    metadata_df <- metadata_df[which(metadata_df[,"Cell_Line"] %in% colnames(PRISM_df)),]

    PRISM_df <- unique(PRISM_df[intersect(colnames(PRISM_df),metadata_df[,"Cell_Line"])])

    
  ##  print(nrow(Achilles_df))
   

  ##  print(nrow(Achilles_df))

    groups <- sort(unique(metadata_df[,metadata_column]))
    grouped_df <- as.data.frame(do.call("cbind", lapply(groups, function(x) rowMeans(PRISM_df[filter(metadata_df, get(metadata_column)==x)[,"Cell_Line"]],na.rm=TRUE))))

                                                                                                                          rownames(grouped_df) <- rownames(PRISM_df)
    colnames(grouped_df) <- groups
    grouped_df <- zscore_df(grouped_df, "row","mean")
    grouped_list <- lapply(colnames(grouped_df), function(x) data.frame(grouped_df[,x],rownames(grouped_df),1-get_quantile(grouped_df[,x]),stringsAsFactors=F,x))
    grouped_list <- lapply(grouped_list, function(x) sort_df(x,1,sort_order=F))
    
    grouped_list <- do.call("rbind", grouped_list)
    colnames(grouped_list) <- c("Zscore","Drug","Quantile","Cancer_Type")
    grouped_list <- split(grouped_list,grouped_list$Cancer_Type)
    grouped_list <- lapply(grouped_list, function(x) filter(x, Quantile>=quantile_cutoff))

##    str(grouped_list)
    

    print("Getting drug targets")

    for(i in names(grouped_list)){ sub <- grouped_list[[i]]
                                   sub_int <- do.call("rbind",lapply(1:nrow(sub), function(x) data.frame(sub[x,],unique(unlist(strsplit(filter(compound_metadata_df, drug_name==sub[x,"Drug"])$target,", "))),stringsAsFactors=F)))
                                   colnames(sub_int)[5] <- "Target"
                                   grouped_list[[i]] <- sub_int
                                   print(paste0("Done getting targets for ",i))

                               }
    grouped_list <- lapply(grouped_list, function(x) sort_df(x,"Zscore",sort_order=FALSE))
    ##    grouped_list <- grouped_list[intersect(unique(carnets_df$Cancer_Type), names(grouped_list))]
    

    print("Running comparisons")

    str(grouped_list)
    
    fisher_list <- lapply(names(grouped_list), function(x) make_and_run_fishers_exact_vectors(unique(filter(carnets_df, Cancer_Type==x)$Gene2), unique(compound_metadata_df$target),unique(grouped_list[[x]][,"Target"]),x,alternative_opt=alternative_opt))
    names(fisher_list) <- names(grouped_list)
    out <- do.call("rbind", fisher_list)
    rownames(out) <- 1:nrow(out)
    out$Signif <- out$P_value <=0.2
    out$P_value <- round(out$P_value, digits=3)

    
    return(out)
                          
    

}

codependency_vs_networks <- function(Achilles_df,TF_target_matrix_dir,metadata_df,metadata_col="Prediction",TF_score_df,subset,genes_of_interest,pval_cutoff=0.2,cor_method=c("spearman","pearson"),precomputed_mat=NULL){

    require(Hmisc)
    require(dplyr)
    require(reshape2)
    print(dim(Achilles_df))
    genes_of_interest <- intersect(genes_of_interest, rownames(Achilles_df))
   
    if(is.null(metadata_df) && is.null(subset)){
        print("Working with all data")
        Achilles_df <- Achilles_df
        flag <- "All"

    } else {
            Achilles_df <- Achilles_df[,intersect(colnames(Achilles_df),metadata_df[,1])]
            metadata_df <- metadata_df[which(metadata_df[,1] %in% colnames(Achilles_df)),]
            index <- unique(unlist(apply(metadata_df,2, function(x) grep(subset,x))))
            print(subset)
            print(index)
            print(metadata_df[index,1])

            ##stop()

            if(length(index) <5 & length(index)>3){ print(paste0("Not enough observations for very reliable analysis, proceeding with subset + ",5-length(index)," random samples"))
                                                    index <- c(index, sample(setdiff(1:nrow(metadata_df), index),5-length(index)))
                                                    flag <- "Subset"
                                                    Achilles_df <- Achilles_df[metadata_df[index,1]]

                                                } else if (length(index) <= 3){ print("0 to 3 cell lines match your search criteria, Running in global correlation mode")
                                                                                ##stop()
                                                                                Achilles_df <- Achilles_df
                                                                                flag <- "Subset"
                                                                } else{
                                                                    Achilles_df <- Achilles_df[metadata_df[index,1]]
    flag <- "Subset"
                                                                }}
    if(is.null(precomputed_mat)==TRUE){
        print("Calculating correlation matrix")
        cormat <- rcorr(t(Achilles_df),type=cor_method)} else if(is.null(precomputed_mat)==FALSE){ print("Working with precomputed matrix")
                                                cormat <- precomputed_mat}
    print("Finished Correlation matrix")

    gc()
    
    cor_vals <- reshape2::melt(cormat[[1]][genes_of_interest,])
    cor_vals[,1:2] <- apply(cor_vals[,1:2],2,as.character)
    print("Finished melting correlation matrix")

    cor_pvals <- reshape2::melt(cormat[[3]][genes_of_interest,])
    cor_pvals[,1:2] <- apply(cor_pvals[,1:2],2,as.character)
    print("Finished melting p-value matrix")
#        str(cor_vals)
#        str(cor_pvals)

    cor_vals$pvals <- cor_pvals[,3]

   sub_vals <- filter(cor_vals,pvals <=pval_cutoff)
    sub_vals$Pair <- paste0(sub_vals[,1],"~",sub_vals[,2])
    str(sub_vals)


    network_df <- data.frame(stringsAsFactors=F)
    TF_score_df <- filter(TF_score_df, Disease==subset)
    drivers <- unlist(strsplit(unique(TF_score_df$Signif_drivers),"-"))
    samples <- unique(TF_score_df$Sample)
    str(samples)


    for(i in genes_of_interest){
        check <- grep(paste0(i,"_"), list.files(TF_target_matrix_dir))
        if(length(check) !=0){
            print(paste0("Working ",i)) 
        sub <- readRDS(paste0(TF_target_matrix_dir,i,"_target_matrix.rds"))
        sub <- sub[samples,]
        sub_summary <- as.data.frame(colSums(sub)/length(samples))
        sub_summary$Gene <- rownames(sub_summary)
            sub_summary$TF <- i
            sub_summary$TF_score_rank <- grep(paste0(i,"$"),drivers)
        colnames(sub_summary)[1] <- "Frequency"
        sub_summary$Pair <- paste0(sub_summary$TF,"~",sub_summary$Gene)
        network_df <- rbind(network_df, sub_summary,stringsAsFactors=F)
        }                         else{ print(paste0("Skipping ",i)) }}

##    str(network_df)


    out <- merge(sub_vals,network_df, by="Pair")
    str(out)


    colnames(out) <- c("Pair","Gene1","Gene2","Correlation","pvalue","Frequency","Gene","TF","TF_score_rank")
    out$Bool <- out$Frequency >= length(samples)/2
    out$Num_Samples <- length(samples)
    out$Cancer_Type <- subset
    out$Method <- cor_method
    out <- out %>% group_by(Bool,TF) %>% mutate(Median_Codependency=median(abs(Correlation),na.rm=TRUE)) %>% data.frame
    
    return(out)
    
}

plot_codep_vs_networks <- function(codep_vs_networks){

    require(dplyr)
    ##codep_vs_networks <- codep_vs_networks[sample(1:nrow(codep_vs_networks),1000),]
    codep_vs_networks <- codep_vs_networks %>% group_by(Bool,Cancer_Type,TF) %>% mutate(Median_Codependency=median(abs(Correlation),na.rm=TRUE)) %>% data.frame

    
    ##hline <- median(abs(codep_vs_networks$Correlation),na.rm=TRUE)
    summary_df <- unique(codep_vs_networks[c("TF","Cancer_Type","Bool","Median_Codependency")])
    summary_df$color <- "gray50"
  ##  str(summary_df)
    summary_out <- data.frame(stringsAsFactors=F)
    for(i in unique(summary_df$TF)) { sub2 <- filter(summary_df,TF==i)

                                      if(abs(filter(sub2, Bool==TRUE)$Median_Codependency) > abs(filter(sub2,Bool!=TRUE)$Median_Codependency)) {  sub2[,"color"] <- "red" } else{ sub2[,"color"] <- "green"}
                                     ## str(sub2)
                                      summary_out <- rbind(summary_out, sub2)
                                                                        
                                  }
    summary_df <- summary_out
##  print(table(summary_df$color))



    p <- ggplot(codep_vs_networks, aes(Bool, abs(Correlation),fill=Bool))+theme_bw()+geom_violin(trim=TRUE,width=0.3)+geom_boxplot(width=0.2, fill="white",color="black")+geom_line(data=summary_df,aes(Bool,abs(Median_Codependency), alpha=0.3,group=TF),color=summary_df$color)+geom_point(data=summary_df,aes(Bool,abs(Median_Codependency),alpha=0.3, group=TF))+facet_wrap(TF~.)
    ##+geom_hline(yintercept=hline, color="black", lwd=1, lty=2, alpha=0.25)

    return(p)}

make_consensus_networks <- function(network_dir,metadata_df,metadata_col, subset,min_fraction=0.5){

    library(igraph)
    samples <- metadata_df[which(metadata_df[,metadata_col] == subset),"Sample"]

    file_list <- list.files(network_dir,".rds",full.names=TRUE)
    file_list <- unique(unlist(lapply(samples, function(x) grep(x, file_list,value=T))))

    mat <- data.frame(0)
    rownames(mat) <- "Empty"
    for(i in file_list) {

        sub <- readRDS(i)
        el <- get.edgelist(sub)
        el <- unique(data.frame(1,paste0(el[,1],"~",el[,2])))
        str(el)
        rownames(el) <- el[,2]
        el[,2] <- NULL
        
        mat <- vector_smartmerge(mat,el)

        print(paste0("Finished ", grep(i, file_list)," of ",length(file_list)))}

    mat <- mat[-1,-1]
    el_out <- rowSums(mat)/ncol(mat)
    el_out <- el_out[which(el_out >=min_fraction)]
    split1 <- unlist(strsplit(names(el_out),"~"))
    TFs <- split1[seq(1,length(split1),2)]
    Genes <- split1[seq(2,length(split1),2)]
    out <- as.matrix(data.frame(TFs,Genes))
    str(out)
    out <- graph_from_edgelist(out,directed=TRUE)


    
    saveRDS(out, paste0(subset,"_consensus_graph.rds"))

    print("Finished")
}

plot_clustering_comparisons <- function(df1, df2, label1,label2,dimensions=c("Dim1","Dim2"),merge_col="Sample",reference_df=c(1,2),cluster_col1="Disease",cluster_col2="Disease",color_scale="Corces_Hex_Codes_Cancer.txt",scale_column="Subtype_ATAC",legend_ncol=2){

    require(dplyr)

    combo_df <- merge(df1,df2,by=merge_col,all=TRUE)
    colnames(combo_df) <- gsub("\\.x",paste0("_",label1),colnames(combo_df))
    colnames(combo_df) <- gsub("\\.y",paste0("_",label2),colnames(combo_df))

    ## str(combo_df)

    if(reference_df==1){ reference_df <- df1
                         ref_col <- cluster_col1
                         ref_label <- label1
                         alt_df <- df2
                         alt_col <- cluster_col2
                         alt_label <- label2
                     } else{
                         reference_df <- df2
                         ref_col <- cluster_col2
                         ref_label <- label2
                         alt_df <- df1
                         alt_col <- cluster_col1
                         alt_label <- label1

                     }

  ##  str(reference_df)
    print(ref_col)
  ##  print(alt_col)
  ##  print(ref_label)
  ##  print(alt_label)
    disease1 <- unique(reference_df[,ref_col])
    disease2 <- unique(alt_df[,alt_col])
  ##  str(disease1)
  ##  str(disease2)

  ##  str(dimensions)


    centroid_df <- as.data.frame(do.call("rbind",lapply(disease1, function(x) colMeans(filter(reference_df, get(ref_col)== x)[,dimensions]))),stringsAsFactors=F)
    centroid_df[,ref_col] <- unique(reference_df[,ref_col])
    ##print(centroid_df)

    centroid_df <- centroid_df[hclust(dist(centroid_df[,dimensions]))$order,]
    factor_levels <- centroid_df[,cluster_col1]
    ref_col_final <- paste0(ref_col,"_",ref_label)
    alt_col_final <- paste0(alt_col,"_",alt_label)
  ##  str(combo_df)
    ##    combo_df[ref_col_final] <- factor(df$Prediction, levels=factor_levels)
##    ref_col_final <- ref_col
##    alt_col_final <- alt_col
    combo_df[ref_col_final] <- factor(combo_df[,ref_col_final], levels=factor_levels)
 

##    print(centroid_df)
##    print(factor_levels)
  

###Quick hack of factor levels for better plots #####

    df_mat <- matrix(0, nrow=length(disease1), ncol=length(disease2))
##    print(dim(df_mat))
    rownames(df_mat) <- disease1
    colnames(df_mat) <- disease2

    summary_df <- as.data.frame(table(combo_df[c(ref_col_final,alt_col_final)]))
    summary_df <- summary_df %>% group_by(get(ref_col_final)) %>% mutate(Freq_Norm=Freq/sum(Freq)) %>% data.frame
    summary_df[,4] <- NULL
##    print(df_mat[1:5,1:5])
##    str(summary_df)
    index <- which(summary_df[,4] !=0) 
    for(n in index){
        i=as.character(summary_df[n,1])
        j=as.character(summary_df[n,2])
        val=summary_df[n,4]
        df_mat[i,j] <- val
    }
##    print(df_mat[1:5,1:5])

    df_mat <- df_mat[as.character(factor_levels),]
    

    j_blacklist <- c()
    indices <- c()
    for(i in 1:nrow(df_mat)){
        best <- sort(df_mat[i,],sort_order=TRUE)
        best_available <- best[setdiff(names(best),j_blacklist)]
        if(best_available[1] >=0.2) {j_blacklist <- c(j_blacklist,names(best_available[1])); indices <- c(indices,rownames(df_mat)[i])} else{ j_blacklist <- j_blacklist; indices <- indices}}
##    print(indices)

    missing <- setdiff(unique(combo_df[,alt_col_final]),j_blacklist)
##    print(missing)
    for(i in missing){
##        print(i)
        best <- sort(df_mat[,i],sort_order=TRUE)
##        print(best)
        if(best[1] !=0){ insert_loc <- grep(names(best)[1], indices)

                         if(length(insert_loc)==0){insert_loc <- grep(names(best)[2], indices)} else{print("Moving to secondary")}

                         str(insert_loc)
                         j_blacklist <- c(j_blacklist[1:insert_loc], i, j_blacklist[(insert_loc+1):length(j_blacklist)])} else{ j_blacklist <- c(j_blacklist,i)}
 ##       print(j_blacklist)

    }


    factor_levels2 <- j_blacklist
    print(factor_levels2)
    combo_df[alt_col_final] <- factor(combo_df[,alt_col_final], levels=factor_levels2)
##    str(combo_df)
   

    p <- ggplot(combo_df, aes_string(alt_col_final,ref_col_final,fill=scale_column))+geom_jitter(shape=21,size=3, width=0.2, height=0.2,show.legend=T)+theme_bw()+theme(axis.text.x=element_text(face="bold",size=8,angle=30,hjust=1),axis.text.y=element_text(face="bold",size=8),strip.text = element_text(colour = "black", face = "bold",size=10))+guides(fill=guide_legend(ncol=legend_ncol))+scale_TF_score(color_list=color_scale,scales="fill")
    print(p)
    print("Finished")

}
 
saturation_curve_vec <- function(list_in){
    ##must be a named list object
    
    list_len <- lapply(list_in, function(x) length(x))

    list_len <- sort(unlist(list_len),decreasing=TRUE)


    initial <- names(list_len)[1]
##    str(initial)
    blacklist <- c(initial)
##    str(blacklist)
    vec_len <- list_in[[initial]]
##    str(vec_len)
    vec_all_len <- length(vec_len)
##    str(vec_all_len)

    for(j in 2:length(list_in)){
        int_lengths <- lapply(list_in, function(x) length(setdiff(x,vec_len)))
##        str(int_lengths)
        sub_length <- sort(unlist(int_lengths[setdiff(names(int_lengths),blacklist)]),decreasing=TRUE)
        index <- names(sub_length)[1]

        vec_len <- union(vec_len, list_in[[index]])
        vec_all_len <- c(vec_all_len, length(vec_len))

        blacklist <- c(blacklist,index)
        
        }
    out <- data.frame(vec_all_len,1:length(vec_all_len),stringsAsFactors=F)
    
    colnames(out) <- c("Union","Num_Samples")
    out$Index <- blacklist
    return(out)

}

plot_IMPACT_mutations <- function(IMPACT_df, subset=NULL,subset_col="TCGA"){
    if(is.null(subset)){
        print("Plotting Mutations for all ")}
    else if (subset =="All") {
        print("Plotting Mutations for all ")}
    else if (subset %in% IMPACT_df[,"TCGA"]){
        print(paste0("Plotting Mutations for ",subset))}
    else{ print(paste0(subset," is not in the provided dataset and column"))}

    


}
    
plot_IMPACT_mut_heatmap <- function(IMPACT_df, mut_col="TCGA_Element_mut_freq2",scale=c(TRUE,FALSE)){
    require(dplyr)
    require(ggplot2)
    require(reshape2)

    blacklist_cancers <- c("SARC","DLBCL","PAAD","Union","OV")
    
    IMPACT_df <- filter(IMPACT_df, Element!="None",!TCGA %in% blacklist_cancers)
##    str(IMPACT_df)
    sub <- unique(IMPACT_df[c(mut_col,"TCGA","Element")])
##    str(sub)

    mat <- matrix(0, nrow=length(unique(sub[,"TCGA"])),ncol=length(unique(sub[,"Element"])))
    rownames(mat) <- unique(sub[,"TCGA"])
    colnames(mat) <- unique(sub[,"Element"])

    for(i in 1:nrow(sub)){ mat[sub[i,2],sub[i,3]] <- sub[i,1]}
    print(mat[1:5,1:5])
    if(scale==TRUE){
        print("Scaling within cancer type")
        mat <- t(apply(mat,1, function(x) x/max(abs(x))))
        dist_mat <- hclust(dist(mat))
    mat <- mat[dist_mat$order,]
    } else {print("Not scaling")
            dist_mat <- hclust(dist(mat))
    mat <- mat[dist_mat$order,]
}
    print(mat[1:5,1:5])

    j_blacklist <- c()
    indices <- c()
    str(dim(mat))

    saveRDS(mat,"Mutation_mat.rds")
    for(i in 1:nrow(mat)){
        best <- sort(mat[i,],decreasing=TRUE)
        best_available <- best[setdiff(names(best),j_blacklist)]
        if(best_available[1] >=0.1) {j_blacklist <- c(j_blacklist,names(best_available[1])); indices <- c(indices,rownames(mat)[i])} else{ j_blacklist <- j_blacklist; indices <- indices}}
    ##    print(indices)
    print(j_blacklist)

    missing <- setdiff(unique(sub[,3]),j_blacklist)
    str(missing)
    for(i in missing){
##        print(i)
        best <- sort(mat[,i],decreasing=TRUE)
##        print(best)
        if(best[1] !=0){ insert_loc <- grep(names(best)[1], indices)

                         if(length(insert_loc)==0){insert_loc <- grep(names(best)[2], indices)} else{print("Moving to secondary")}

                         str(insert_loc)
                         j_blacklist <- c(j_blacklist[1:insert_loc], i, j_blacklist[(insert_loc+1):length(j_blacklist)])} else{ j_blacklist <- c(j_blacklist,i)}
 ##       print(j_blacklist)

    }


    

##    factor_levels <- centroid_df[,cluster_col1]



}


seriate_matrix <- function(mat, reference_margin=c("row","column"),min_cutoff=0.1){

    sub <- mat
    j_blacklist <- c()
    indices <- c()

    if(reference_margin=="row"){
        alt_list <- colnames(mat)
        index <- hclust(dist(sub))$order
##        print(rownames(sub)[index])
        sub <- sub[index,]
    for(i in 1:nrow(sub)){
        best <- sort(sub[i,],decreasing=TRUE)
        best_available <- best[setdiff(names(best),j_blacklist)]
        if(best_available[1] >=min_cutoff) {j_blacklist <- c(j_blacklist,names(best_available[1])); indices <- c(indices,rownames(sub)[i])} else{ j_blacklist <- j_blacklist; indices <- indices}}
##       print(indices)
##    print(j_blacklist)

    missing <- setdiff(alt_list,j_blacklist)
##    str(missing)
    for(i in missing){
    ##    print(i)
        best <- sort(sub[,i],decreasing=TRUE)
        ##print(best)
        if(best[1] !=0){ insert_loc <- grep(names(best)[1], indices)
  ##                       print("Primary:")
       
                         if(length(insert_loc)==0){insert_loc <- grep(names(best)[2], indices); print("Moving to secondary") } 
        
        j_blacklist <- c(j_blacklist[1:insert_loc], i, j_blacklist[(insert_loc+1):length(j_blacklist)])} else{ j_blacklist <- c(j_blacklist,i); print("Can't find a position for this element")}}
        sub <- sub[,j_blacklist]

    } else if (reference_margin=="column") {

        alt_list <- rownames(sub)
        index <- hclust(dist(t(sub)))$order
        ##        print(colnames(sub)[index])

        sub <- sub[,index]

            for(i in 1:ncol(sub)){
                best <- sort(sub[,i],decreasing=TRUE)

        best_available <- best[setdiff(names(best),j_blacklist)]
##        str(best_available)
##        str(length(best_available))
        if(best_available[1] >=min_cutoff && length(best_available)!=0) {j_blacklist <- c(j_blacklist,names(best_available[1])); indices <- c(indices,colnames(sub)[i])} else { print("No value found"); j_blacklist <- j_blacklist; indices <- indices}; }
##        print(indices)
##    print(j_blacklist)

        missing <- setdiff(alt_list,j_blacklist)
##        print("Missing:")
    str(missing)
        for(i in missing){
  ##          print(i)
            best <- sort(sub[i,],decreasing=TRUE)
  ##          str(best)
        if(best[1] !=0){ insert_loc <- grep(names(best)[1], indices)

                         if(length(insert_loc)==0){insert_loc <- grep(names(best)[2], indices)} else{print("Moving to secondary")}

  ##                       str(insert_loc)
                         j_blacklist <- c(j_blacklist[1:insert_loc], i, j_blacklist[(insert_loc+1):length(j_blacklist)])} else{ j_blacklist <- c(j_blacklist,i)}}
        sub <- sub[j_blacklist,]

    }


    return(sub)
        
}


    
is.directory <- function(obj){
    status <- file.info(obj)$isdir
    return(status)}




set_walltime <- function(.Object,time=NULL)
{
    if(is.null(time)){
        bcmds <- .Object@runinfo$bcmd
        cmd_check <- grep("[0-9][0-9]:[0-9][0-9]", bcmds[1])
        if(length(cmd_check) == 0){
            bcmds <- unlist(lapply(bcmds, function(x) gsub("  -R","  -W 12:00 -R",x)))} else{
                                                                                          bcmds <- unlist(lapply(bcmds, function(x) gsub("[0-9][0-9]:[0-9][0-9]","12:00",x)))}
        
        .Object@runinfo$bcmd <- bcmds
    } else{
                       time_check <- grep("^[0-9][0-9]:[0-9][0-9]$|^[0-9][0-9][0-9]:[0-9][0-9]$", time)
                       
                                                                                      if(length(time_check) ==1){
                                                                                          bcmds <- .Object@runinfo$bcmd
                                                                                          cmd_check <- grep("[0-9][0-9]:[0-9][0-9]", bcmds[1])
                                                                                                  if(length(cmd_check) == 0){
                                                                                          bcmds <- unlist(lapply(bcmds, function(x) gsub("  -R",paste0("  -W ",time," -R"),x)))} else{
                                                                                          bcmds <- unlist(lapply(bcmds, function(x) gsub("[0-9][0-9][0-9]:[0-9][0-9]|[0-9][0-9]:[0-9][0-9]",time,x)))}
                                                                                          .Object@runinfo$bcmd <- bcmds
                                                                                      } else{ print("Malformed walltime. Not updating"); stop()}}
return(.Object)
}


        
get_walltime <- function(.Object){
    require(stringr)
    bcmds <- .Object@runinfo$bcmd
    times <- unlist(lapply(bcmds, function(x) str_extract(x,"[0-9][0-9]:[0-9][0-9]")))
    names(times) <- names(bcmds)
    print(times)}



alignment_qc <-  function(indir,sample_ID){
    require(Rsamtools)
    
    infiles <- list.files(indir)

    sam_files <- infiles[grep("*.sam",infiles)]
    bam_files <- infiles[grep("*.bam$",infiles)]
    trimmed_fastq <- infiles[grep("_trimmed.fastq$",infiles)]
    fastq <- setdiff(infiles[grep(".fastq$", infiles)], trimmed_fastq)

    if(!all(unlist(lapply(list(sam_files, bam_files,trimmed_fastq,fastq), function(x) length(x)>=1)))){ print("Missing one or more filetypes"); stop()} else { print("All filetypes present and accounted for")}


    system("set -euxo pipefail")
    ##print(sam_info)
    bam_info <- try(do.call("rbind",lapply(bam_files, function(x) data.frame(x, file.size(paste0(indir,x))/1e9, as.numeric(system(paste0("samtools view -c ",indir,x),intern=TRUE))/2,"BAM"))))
    ##print(bam_info)
    print("Finished bam")
    sam_info <- try(do.call("rbind",lapply(sam_files, function(x) data.frame(x, file.size(paste0(indir,x))/1e9, as.numeric(system(paste0("samtools view -c ",indir,x),intern=TRUE))/2,"SAM"))))
    print("Finished sam")

    trimmed_fastq_info <- try(do.call("rbind",lapply(trimmed_fastq, function(x) data.frame(x, file.size(paste0(indir,x))/1e9, as.numeric(unlist(strsplit(system(paste0("wc -l ", indir,x), intern=TRUE)," "))[1])/4,"TRIMMED_FQ"))))
    print("Finished trimmed_fastq")
    fastq_info <- try(do.call("rbind",lapply(fastq, function(x) data.frame(x, file.size(paste0(indir,x))/1e9, as.numeric(unlist(strsplit(system(paste0("wc -l ", indir,x), intern=TRUE)," "))[1])/4,"FQ"))))
    print("Finished source fastq")
    


    colnames(sam_info) <- colnames(bam_info) <- colnames(trimmed_fastq_info) <- colnames(fastq_info) <- c("File","Size","Reads","Type")
    
    out <- do.call("rbind",list(sam_info,bam_info, trimmed_fastq_info, fastq_info))
    out$Sample <- sample_ID
    ##print(out)
##    df <- as.data.frame(table(scanBam(paste0(indir,bam_files))[[1]]$mapq))
##    out_list <- list(out, df)
    ##    return(out_list)
    out_list <- list(out)
    return(out_list)
}

##    
##    str(file.size(paste0(indir,fastq)))
##    
##    names(out_list) <- c("SAM","BAM","TRIMMED_FQ","FQ")

##    for(i in 1:length(out_list)){
##        sub <- out_list[[i]]
##        if(length(sub) !=0) { df <- do.call("rbind", lapply(sub, function(x) data.frame(x, file.size(paste0(indir,x)), as.numeric(unlist(strsplit(system(paste0("wc -l ", indir,x), intern=TRUE)," "))[1]))))
##            print(df)
            


parse_featureCounts <- function(featureCount_df,drop_haplotypes=TRUE){
    if(is.character(featureCount_df)){
        featureCount_df <- read.table(featureCount_df, sep='\t', header=T, stringsAsFactors=F)} else if(class(featureCount_df)=="data.frame"){
                                                                                                  featureCount_df <- featureCount_df}
    out <- list()

    for(i in 1:nrow(featureCount_df)){
##    for(i in 1:100){
        sub <- featureCount_df[i,]
        int <- as.data.frame(do.call("cbind",lapply(sub, function(x) unlist(strsplit(as.character(x),";")))))

        out <- c(out, list(int))

}

    out <- do.call("rbind",out)


##    str(out)
    if(drop_haplotypes==TRUE){
        index1 <- grep("^chr[0-9]{1,2}$",out$Chr)
        index2 <- grep("^chr[A-Z]{1,2}$", out$Chr)
        
        out <- out[c(index1,index2),]
    }
    out <- unique(out[,c(2,3,4,5,1,7)])
    rownames(out) <- NULL
    colnames(out)[6] <- "Count"

    return(out)
    }
make_featureCount_matrix <- function(dir, extension="_featureCounts$",file_list=NULL){
    source("~/Andre_F_functions.R")
    require(dplyr)
    if(is.null(file_list)){
        file_list <- list.files(dir, extension)}
    else{ file_list <- file_list}
        IDs <- gsub(extension,"", file_list)
        out <- parse_featureCounts(paste0(dir,file_list[1]))
        colnames(out)[6] <- IDs[1]
        out[,6] <- as.numeric(out[,6])
        for(i in 2:length(file_list)){
            int <- parse_featureCounts(paste0(dir,file_list[i]))
            out[IDs[i]] <- as.numeric(int$Count)
            print(paste0("Finished ",i," of ", length(file_list)))}
        return(out)

    }
make_nucleo_job <- function(sample_manifest, working_dir,nuc_json="/work/bergerm1/Innovation/share/nucleo/nucleo_blank_job.json", nuc_sh="/work/bergerm1/Innovation/share/nucleo/nucleo_blank_run_cmd.sh",ignore_missing=FALSE,outdir="~/BergerLab_Work/") {
    require(dplyr)
    ##sample manifest needs to be a headerless tsv of format sample_id \t fastq1 \t fastq2, unique sample_ids are required for each line and full paths for fastqs

    sample_manifest <- read.table(sample_manifest, header=F, sep='\t',stringsAsFactors=F)
    sample_manifest <- unique(sample_manifest)
    manifest_format_check <- ncol(sample_manifest)
    if (manifest_format_check ==3){
        manifest_check <- unlist(lapply(1:nrow(sample_manifest), function(x) paste0(file.exists(sample_manifest[x,2]),"~",file.exists(sample_manifest[x,3]))))
        sample_manifest$check <- manifest_check

        if(any(sample_manifest$check != "TRUE~TRUE")){ print("Some provided files in manifest do not exist. Missing file(s):")
            print(filter(sample_manifest, check!="TRUE~TRUE"))
            if(ignore_missing==FALSE){ stop()} else{ print("Removing samples with missing files and proceeding.")
                                                        sample_manifest <- filter(sample_manifest, check=="TRUE~TRUE")}

        } else{ print("All files in manifest accounted for."); sample_manifest <- sample_manifest}
        }
    else{ print("Malformed manifest file"); stop()}

    rownames(sample_manifest) <- sample_manifest[,1]


    for(i in unique(sample_manifest[,1])){
        filled_json <- gsub("sample_name", i, readLines(nuc_json))
        filled_json <- gsub("path/to/fq1", sample_manifest[i,2],filled_json)
        filled_json <- gsub("path/to/fq2", sample_manifest[i,3],filled_json)
        writeLines(filled_json, paste0(outdir,i,"_job.json"))}

    out_script <- c()
    for(i in unique(sample_manifest[,1])){
        print(i)
        script <- readLines(nuc_sh)
        str(script)
        filled_script <- gsub("sample_name",i, script)
        filled_script <- gsub("path/to/wkdir", working_dir,filled_script)
        out_script <- c(out_script,filled_script)
        }
    writeLines(out_script,paste0(outdir,"all_run_nucleo.sh"))
str(nrow(sample_manifest))
    print(paste0("Finished prepping nucleo jobs for ",nrow(sample_manifest)," samples"))

    }

CCS_fusions_to_bedpe <- function(CCS_fusions,add_chr=TRUE,name_root="Fusion_"){
    if(is.character(CCS_fusions)){
        CCS_fusions <- read.delim(CCS_fusions,header=T, sep='\t', stringsAsFactors=F)} else if (class(CCS_fusions) == "data.frame"){ CCS_fusions <- CCS_fusions} else{ print(" Check input fusions. Not a dataframe or parsable table")}

##    str(CCS_fusions)
    region1 <- do.call("rbind",lapply(CCS_fusions$breakpoint1, function(x) as.data.frame(cbind(unlist(strsplit(x,":"))[1],as.numeric(unlist(strsplit(x,":"))[2]),as.numeric(unlist(strsplit(x,":")))[2]+1))))
    region2 <- do.call("rbind",lapply(CCS_fusions$breakpoint2, function(x) as.data.frame(cbind(unlist(strsplit(x,":"))[1],as.numeric(unlist(strsplit(x,":"))[2]),as.numeric(unlist(strsplit(x,":")))[2]+1))))

    if(add_chr ==TRUE){
        region1[,1] <- paste0("chr",region1[,1])
        region2[,1] <- paste0("chr",region2[,1])
    }
    colnames(region1) <- c("chr_bp1","start_bp1","end_bp1")
    colnames(region2) <- c("chr_bp2","start_bp2","end_bp2")



    CCS_fusions <- cbind(region1, region2, CCS_fusions)

    CCS_fusions <- bedpe_to_granges(CCS_fusions,el_names=name_root)
   ## print(CCS_fusions)
    
    return(CCS_fusions)
    
    }

make_consensus_peakset <- function(dir,extension,subdirs=TRUE,filter_col=NULL,filter_cutoff=NULL,fix_width=NULL,file_list=NULL){
    require(data.table)
    require(gUtils)
    require(GenomicRanges)


    if(!is.null(file_list)){
        file_list <- file_list} else{
    if(subdirs ==TRUE){

       file_list <- system(paste0("find ",dir," -name \'*",extension,"\' "),intern=TRUE)
##        str(file_list)
    } else {
        file_list <- list.files(dir, extension,full.names=TRUE)
##        str(file_list )
    }}

    file_list <- sample(file_list)
##   file_list <- file_list[1:5]
    ID_list <- unlist(lapply(file_list, function(x) gsub(extension,"",data.table::last(unlist(strsplit(x,"/"))))))
    str(ID_list)
    str(file_list)
    print(file_list[1:5])
    grl <- lapply(file_list, function(X) bed_to_granges_dynamic(X,header=FALSE))
##    print(grl[[1]])
    if(!is.null(filter_col)){
        grl <- lapply(grl, function(x) x[which(mcols(x)[,filter_col] >= filter_cutoff)])} else { grl <- grl}

    if(!is.null(fix_width)){
        grl <- lapply(grl,(function(x) gr.mid(x)+(fix_width/2)))} else{ grl <- grl}
    
##        print(grl[[1]])
    names(grl) <- ID_list

  
    peak_set <- unique(grl[[1]])
    peak_set$count <- 1

    for(i in 2:length(grl)){
        print(paste0("Working on ", names(grl)[i]))
        missing <- setdiff_with_metadata(grl[[i]],peak_set)

        similar <- findOverlaps(peak_set,grl[[i]])@from
##        print(length(similar))

        if(length(similar) >=1){ peak_set[similar]$count <- peak_set[similar]$count+1} 
        

        if(length(missing) >=1){
            missing$count <- 1
            peak_set <- unique(c(peak_set, missing))} else { peak_set <- peak_set}
        print(summary(peak_set$count))
        print(paste0("Finished ",i," of ", length(grl)))

    }
    

    return(peak_set)
    

    }

generate_peak_height_matrix <- function(dir, extension, subdirs=TRUE,regions,outfile,file_list=NULL,compressed_output="multibig.npz"){
    if(!is.null(file_list)){
        file_list <- file_list
        str(file_list)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
        str(labels)
    } else{

    if(subdirs ==TRUE){
        file_list <- system(paste0("find ",dir," -name \'*",extension,"\' "),intern=TRUE)
        str(file_list)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
        str(labels)
    } else {
        file_list <- list.files(dir, extension,full.names=TRUE)
        str(file_list)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
        str(labels)
        
    }}

##    file_list <- file_list[1:5]
    if(is.character(regions)){
        cmd <- paste0("multiBigwigSummary BED-file -b ", paste(file_list, collapse=" ")," -o ", compressed_output," --outRawCounts ", outfile," --BED ", regions," --labels ",paste(labels,collapse=" "))} else if(class(regions) =="GRanges") {
                                                                                                                                                                                                  gr_to_bed(regions,"func.bed",metadata=TRUE)
                                                                                                                                                                                                  cmd <- paste0("multiBigwigSummary BED-file -b ", paste(file_list, collapse=" ")," -o ",compressed_output," --outRawCounts ", outfile," --BED func.bed --labels ",paste(labels,collapse=" "))
                                                                                                                                                                                              }
##    print(cmd)
   system(cmd)
}


generate_peak_summit_matrix <- function(dir, extension=".bigwig", subdirs=TRUE, regions, file_list=NULL,labels=NULL,parallel=FALSE, n_cores=5, outfile=NULL,verbose=TRUE,operation=c("max","median","mean")){
    require(rtracklayer)

    
        if(!is.null(file_list)){
        file_list <- file_list

        if(is.null(labels)){
        labels <- gsub(paste0(dir,"|",extension),"",file_list)} else{ labels <- labels}

        } else{
                    if(is.null(labels)){
        labels <- gsub(paste0(dir,"|",extension),"",file_list)} else{ labels <- labels}


    if(subdirs ==TRUE){
        file_list <- system(paste0("find ",dir," -name \'*",extension,"\' "),intern=TRUE)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
    } else {
        file_list <- list.files(dir, extension,full.names=TRUE)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
    }}
##    str(file_list)
    str(labels)
    if(is.character(regions)){
        regions <- bed_to_granges_dynamic(regions)
        names(regions) <- regions$V4
        print("Importing intervals from bed")

        } else if(class(regions) =="GRanges") {
            print("Working with the provided intervals")}
##        print(regions)
##    regions <- regions[1:10]
        grl <- split(regions, names(regions))

    out_mat <- data.frame(rep(0, length(regions)),row.names=names(regions))
    colnames(out_mat) <- "Empty"

    labels <- sort(labels)
    if(parallel==FALSE){
    for(i in labels){
        print(paste0("Working ",i))
        infile <- grep(i, file_list, value=T)
##        print(infile)

        inwig <- import.bw(infile, selection=regions+0.5e3)
        ##        print(length(inwig))
        if(operation=="max"){
                     summits <- unlist(lapply(grl, function(x) max(unique(intersect_with_metadata(inwig,x))$score,na.rm=TRUE)))} else if(operation=="mean"){summits <- unlist(lapply(grl, function(x) mean(unique(intersect_with_metadata(inwig,x))$score,na.rm=TRUE)))} else if(operation=="median"){summits <- unlist(lapply(grl, function(x) median(unique(intersect_with_metadata(inwig,x)$score),na.rm=TRUE)))}
                    


        if(verbose){
        str(summits)}

        summits <- summits[rownames(out_mat)]
        out_mat[,i] <- summits


    }
        out_mat[,1] <- NULL
    } else if(parallel==TRUE){
         print("Initializing parallel")
         require(foreach)
         require(doMC)
         registerDoMC(cores=n_cores)
         source("/home/forbesa1/Andre_F_functions.R")

        out_mat <- foreach(i=labels,.combine="cbind",.multicombine=TRUE) %dopar% {
                     infile <- grep(paste0(i,extension), file_list, value=T)

                     inwig <- import.bw(infile, selection=regions+2.5e3)
                     ## print(inwig)
                     if(operation=="max"){
                     summits <- unlist(lapply(grl, function(x) max(intersect_with_metadata(inwig,x)$score)))} else if(operation=="mean"){summits <- unlist(lapply(grl, function(x) mean(intersect_with_metadata(inwig,x)$score)))} else if(operation=="median"){summits <- unlist(lapply(grl, function(x) median(intersect_with_metadata(inwig,x)$score)))}
                    if(verbose){str(summits)}

                     print(i)
                     

        summits <- summits[rownames(out_mat)]
        return(summits)
                     
        }

         colnames(out_mat) <- labels
         str(out_mat)}
                    


    out_mat <- as.data.frame(out_mat)
    if(is.null(outfile)){
        return(out_mat)} else{ saveRDS(out_mat,outfile)}
        
       }
                                          
generate_peak_binary_matrix <- function(dir, extension, subdirs=TRUE, regions, file_list=NULL,region_names="name",filter_col="V7",filter_val=5){
    require(rtracklayer)

    
        if(!is.null(file_list)){
        file_list <- file_list

        labels <- gsub(paste0(dir,"|",extension),"",file_list)

    } else{

    if(subdirs ==TRUE){
        file_list <- system(paste0("find ",dir," -name \'*",extension,"\' "),intern=TRUE)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
    } else {
        file_list <- list.files(dir, extension,full.names=TRUE)
        labels <- gsub(paste0(dir,"|",extension),"",file_list)
    }}


        if(is.character(regions)){
            regions <- bed_to_granges_dynamic(regions)
            out_mat <- data.frame(rep(0,length(regions)), row.names=mcols(regions)[,region_names])
            names(regions) <- mcols(regions)[,region_names]
##            print(regions)
        print("Importing intervals from bed")

        } else if(class(regions) =="GRanges") {
            print("Working with the provided intervals")
 
##            print(regions)
            out_mat <- data.frame(rep(0,length(regions)), row.names=names(regions))}

    colnames(out_mat) <- "Empty"

    file_list <- sort(file_list)
##    file_list <- file_list[1:11]
    for(i in file_list){
        
        ID <- gsub(paste0(dir,"|",extension),"",i)
        ##print(ID)
        peaks <- bed_to_granges_dynamic(i)
        peaks  <- filter_gr(peaks, filter_col, filter_val, "greater_equal")
        ##print(peaks)
        int <- names(intersect_with_metadata(regions, peaks))
        out_mat[,ID] <- 0
        out_mat[int,ID] <- 1

    }
    out_mat <- out_mat[,-1]
    return(out_mat)

    }



generate_expression_matrix_kallisto <- function(dir, extension="abundance.tsv", subdirs=TRUE,file_list=NULL,abundance_col="tpm"){
    if(!is.null(file_list)){
        file_list <- file_list} else{

    if(subdirs ==TRUE){
        file_list <- system(paste0("find ",dir," -name \'*",extension,"\' "),intern=TRUE)
  ##      str(file_list)
    } else {
        file_list <- list.files(dir, extension,full.names=TRUE)
  ##      str(file_list)
        
    }}



    mat <- data.frame(0)
    rownames(mat) <- "Empty"
    for(i in file_list){
        df <- read.delim(i, sep='\t', stringsAsFactors=F,header=T)
##        df <- df[1:1000,]
        df <- df[c("target_id",abundance_col)]
        rownames(df) <- df[,1]
##        str(df)
        df[,1] <- NULL
        mat <- vector_smartmerge(mat,df)
        
        }

    mat <- mat[-1,-1]
    colnames(mat) <- gsub(paste0("Kallisto/|kallisto/|",extension),"",file_list)

    return(mat)

    }
                                               

gr_remove_chromosomes <- function(gr,valid_list=NULL){
    if(is.null(valid_list)){
        valid_chr <- c(1:22,"X","Y","M","MT")
        valid_chr <- c(valid_chr, paste0("chr",valid_chr))
        } else {valid_chr <- valid_chr }


    gr <- gr[which(seqnames(gr) %in% valid_chr)]

    seqlevels(gr) <- as.character(unique(seqnames(gr)@values))

    return(gr)
}
    
trim_invalid_intervals <- function(gr, genome=c("hg19","hg38"),strand_col=NULL){

    require(dplyr)
    require(reshape2)
    if(genome=="hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19)} else if(genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38)}
    source("~/Andre_F_functions.R")

    lens <- seqlengths(Hsapiens)
    valid_gr <- data.frame(names(lens), 1, lens)

    rownames(valid_gr) <- names(lens)

    gr <- gr_to_bed(gr, outfile=NULL, TRUE)
    gr$Max <- valid_gr[gr$Chrom,"lens"]

##    str(gr)
    index_start1 <- which(gr$Start <0)
    index_start2 <- which(gr$Start > gr$Max)

    index_end <- which(gr$End > gr$Max)

##    str(index_start1)
##    str(index_start2)
##    str(index_end)
    

    gr[index_start1,"Start"] <- 1
    gr[index_start2, "Start"] <- gr[index_start2,"Max"]
    gr[index_end,"End"] <- gr[index_end,"Max"]

    gr$Max <- NULL
    

    if(is.null(strand_col)){
        print("no strand info provided")
    gr <- table_to_granges(gr)} else { print("strand info provided"); table_to_granges(gr, strand_col)}
    return(gr)
    }
subtract_gr <- function(gr1, gr2,collapse=FALSE,stranded=FALSE,keep_metadata=TRUE){
    ##Hacky and requires system call to installed bedtools. Will clobber files if run in parallel
    
    gr_to_bed(gr1,"gr1.bed")
    gr_to_bed(gr2,"gr2.bed")
    if(stranded){system("bedtools subtract -s -a gr1.bed -b gr2.bed > gr3.bed")} else{
    system("bedtools subtract -a gr1.bed -b gr2.bed > gr3.bed")}
    gr_out <- bed_to_granges("gr3.bed")

##    print(gr_out)
   ## gr1 <- reduce(gr1)
   ## print(gr1)
   ## gr2 <- reduce(gr2)
   ## print(gr2)
###z <- disjoin(c(gr1, gr2))
    if(collapse==TRUE) { gr_out <- reduce(gr_out)} else{ gr_out <-gr_out}
    ##    print(gr_out)

    if(keep_metadata){
        print("This feature is a WIP, May be buggy")
##        print(gr_out)
##        print(gr1)
       
        overlaps <- findOverlaps(gr_out,gr1)
##        print(overlaps)
        meta <- mcols(gr1)[overlaps@to,]
        
##        str(meta)
        mcols(gr_out) <- meta
        strand(gr_out) <- strand(gr1)[overlaps@to]}

    end(gr_out) <- end(gr_out)-1L

    lapply(1:3, function(x) system(paste0("rm gr",x,".bed")))
        
    return(gr_out)

}


make_anf_anno_from_gtf <- function(gtf,group_column="type",genome=c("hg19","hg38")){
    require(rtracklayer)
    require(GenomicRanges)
    if(genome=="hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19)} else if(genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38)}
    source("~/Andre_F_functions.R")

    if(is.character(gtf)){ gtf <- rtracklayer::import(gtf)} else if(class(gtf) == "GRanges"){
                                                              gtf <- gtf}
    
    anno_types <- unique(mcols(gtf)[,group_column])
    print(paste0("Using annotation types:", paste(anno_types, collapse=" ")))

    anno_list <- lapply(anno_types, function(x) gtf[which(mcols(gtf)[,group_column] ==x)])
    names(anno_list) <- anno_types
    anno_list$Intron <- subtract_gr(anno_list[[grep("transcript", names(anno_list))]],anno_list[[grep("exon", names(anno_list))]])
    anno_list$Intron$type <- "Intron"
    anno_list$Promoter <- get_promoters(anno_list[[grep("*start*", names(anno_list))]], 1000)
    anno_list$Promoter$type <- "Promoter"
    for(i in 1:length(anno_list)){
        sub <- anno_list[[i]]
        seqlevels(sub) <- seqlevels(Hsapiens)
        seqinfo(sub) <- seqinfo(Hsapiens)
        anno_list[[i]] <- sub }

    anno_list$Intergenic <- unique(subtract_gr(gaps(anno_list[[grep("transcript", names(anno_list))]]),anno_list[[grep("transcript", names(anno_list))]]))
    anno_list$Intergenic$type <- "Intergenic"

##    print(anno_list)

    return(anno_list)
    }


filter_gr <- function(gr,column,column_val,comparison=c("greater","lesser","greater_equal","lesser_equal","equal","in")) {
    if(comparison=="greater"){ out_gr <- gr[which(mcols(gr)[,column] > column_val)] } else if(comparison=="lesser"){ out_gr <- gr[which(mcols(gr)[,column] < column_val)]} else if(comparison=="greater_equal") { out_gr <- gr[which(mcols(gr)[,column] >= column_val)] } else if (comparison=="lesser_equal"){ out_gr <- gr[which(mcols(gr)[,column] <= column_val)]} else if (comparison=="equal"){  out_gr <- gr[which(mcols(gr)[,column]== column_val)]} else if(comparison=="in"){ out_gr <- gr[which(mcols(gr)[,column] %in% column_val)]}

    return(out_gr)}
    

multiway_comparison <- function(peak_list,n_peaks_sampled=2){
    require(gplots)
                                        #Input should be a named list containing 2 or more vectors. The output is a sub-sampled list containing elements present in each comparison of the initial list.
    peak_list <- lapply(peak_list, unique)
    intersections <- attr(venn(peak_list, show.plot=FALSE),"intersections")
##    str(intersections)

    if(!is.null(n_peaks_sampled)){
    out <- lapply(intersections, function(x) sample(x,min(n_peaks_sampled,length(x))))
    }
    else{ out <- intersections}

    return(out)
}


is.nested <- function(x) {
  stopifnot(is.list(x))
  any(sapply(x, function(x) any(class(x) == "list")))
}


run_randomForest <- function(matrix, metadata_df, metadata_column="Disease", sample_margin=c("column","row"),...){

    require(data.table)
    require(randomForest)
    if(class(matrix) == "matrix"){
        if(sample_margin == "row"){ matrix <- as.data.frame(matrix); print("Initializing")} else if(sample_margin =="column"){ matrix <- as.data.frame(t(matrix))
                                                             print("Transposing Matrix then Initializing")}} else if(class(matrix) %in% c("data.table","array","data.frame")){ matrix <- matrix
                                                                                                               if(sample_margin == "row"){ print("Coercing to Matrix then Initializing")} else if(sample_margin =="column"){ matrix <- t(matrix)
                                                                                                                                                                                            print("Coercing to Matrix, Transposing Matrix then Initializing")}}
##    str(matrix)
    
##    str(metadata_df)
    index <- grep("sample", ignore.case=TRUE,colnames(metadata_df))[1]
    metadata_df <- metadata_df[which(metadata_df[,index] %in% rownames(matrix)),]
    
    matrix <- matrix[metadata_df[,index],]
##    str(matrix)
    meta <- metadata_df[rownames(matrix),metadata_column]
    
    matrix <- cbind(matrix,as.factor("Empty"))
    matrix <- as.data.frame(matrix)                                      
    matrix[,ncol(matrix)] <- as.factor(meta)


    colnames(matrix)[ncol(matrix)] <- metadata_column
    print(class(matrix[,ncol(matrix)]))
    str(matrix[,ncol(matrix)])

##    str(matrix)
    rownames(matrix) <- gsub("-","_",rownames(matrix))
        ##str(sub_df)
##        sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),] 
##    transposed_df <- as.data.frame(t(sub_df))
##    transposed_df[,cluster_column] <- as.factor(tsne_df[colnames(expression_df),cluster_column])
##    uniq_classes <- unique(tsne_df[,cluster_column])
##    print(index)
##   print(cluster_column)
##    str(transposed_df)
    print("Starting Random Forest")
##    if(is.null(cutoff)){  cutoff <- 1/length(unique(metadata_df[,metadata_column])); cutoff <- rep(cutoff,length(unique(metadata_df[,metadata_column])))} else { cutoff <- cutoff}
##    print(cutoff)
    rf <- randomForest(x=matrix[,1:(ncol(matrix)-length(metadata_column))],y=matrix[,metadata_column],data=matrix,importance=TRUE,...)
##        importance_df <- as.data.frame(importance(rf))[setdiff(colnames(importance(rf)),uniq_classes)]
##      

    return(rf)



    }
    

bam_read_length_summary <- function(bam, regions=NULL, outfile=ribinNULL,delim="\\n"){

    if(!is.null(regions)){
        if(is.character(regions)){
            minibam_cmd <- paste0("samtools view -b -M -L ", regions," ", bam," > ~/minibam.bam")


        } else if(class(regions) =="GRanges") {
            gr_to_bed(regions,"bam_read_len.bed",metadata=TRUE)
            minibam_cmd <- paste0("samtools view -b -M -L bam_read_len.bed ", bam," > ~/minibam.bam")}

        print("Subsetting bam file to regions of interest")
        system(minibam_cmd)
        print("Indexing shrunken bam file")
        system(paste0("samtools index -b ~/minibam.bam"))
        bam <- "~/minibam.bam"
    }
    


    if(is.null(outfile)){
        print("Getting fragment length distribution")

        cmd <- paste0("bamPEFragmentSize --bamfiles ", bam, " --outRawFragmentLengths ~/frag_len")
                print(cmd)

        system(cmd)
        out <- read.table("~/frag_len", header=T, sep='\t', stringsAsFactors=F)
        return(out)
        
    } else {
        print("Getting fragment length distribution")
        cmd <- paste0("bamPEFragmentSize --bamfiles ", bam, " --outRawFragmentLengths ", outfile)

        system(cmd)}
print("Done")
}

get_zero_coverage_regions_bam <- function(bam,unit=1e9,valid_chr=NULL,focus=c("zero","covered")){
    if(is.null(valid_chr)){
        valid_chr <- c(1:22,"X","Y","M","MT")
        valid_chr <- c(valid_chr, paste0("chr",valid_chr))}
    valid_chr <- paste0("'^(", paste0(valid_chr,collapse="|"),")\\>'")

    out_cov <- gsub(".bam","_zero_cov.txt",bam)
    out_cov2 <- gsub(".bam","_zero_cov2.txt",bam)
    print("1")
    if(focus=="zero"){
        system(paste0("bedtools genomecov -ibam ", bam," -bga | awk '$4==0' > ",out_cov))} else if(focus=="covered"){
                                                                                                     system(paste0("bedtools genomecov -ibam ", bam," -bga | awk '$4!=0' > ",out_cov))}
    system(paste0("grep -E ", valid_chr, "< ", out_cov," > ",out_cov2))
    zero_cov <- as.numeric(system(paste0("awk '{sum += $3 - $2} END {print sum }' ", out_cov2),intern=TRUE))

    system(paste0("rm ",out_cov))
    system(paste0("rm ",out_cov2))
    zero_cov <- zero_cov/unit
    return(zero_cov)}
    

estimate_coverage_bam <- function(bam, regions,subset=TRUE,sample_size=1e4){

    out_bam <- gsub(".bam","_minibam.bam",bam)
    out_cov <- gsub(".bam","_cov.txt",bam)
    if(!is.null(regions)){
        if(is.character(regions)){
            print("Working with bed file")
            region_gr <- reduce(bed_to_granges_dynamic(regions))
            region_gr <- gr_remove_chromosomes(region_gr)

            if(is.null(sample_size)){ region_gr <- region_gr } else{ region_gr <- region_gr[sample(1:length(region_gr), sample_size)]}
            minibam_cmd <- paste0("samtools view -b -M -L ", regions," ", bam," > ", out_bam)


        } else if(class(regions) =="GRanges") {
##            print(regions)
            print("Working with GRanges object")
            region_gr <- reduce(regions)
            
            if(is.null(sample_size)){region_gr <- region_gr } else{ region_gr <- region_gr[sample(1:length(region_gr), sample_size)]}
            random_num <- sample(1:10e4,1)
            gr_to_bed(region_gr,paste0("est_coverage",random_num,".bed"))
            minibam_cmd <- paste0("samtools view -b -M -L est_coverage",random_num,".bed ",bam," > ", out_bam )}
##        print(minibam_cmd)
        system(minibam_cmd)

        system(paste0("bedtools genomecov -ibam ",out_bam," -bga > ",out_cov ))

        cov_file <- bed_to_granges_dynamic(out_cov)
        cov_file <- unique(gr_remove_chromosomes(cov_file))
        ##sub <- reduce(cov_file)
        
        if(subset==TRUE){
            width(cov_file) <- width(cov_file)-1
            
            blacklist <- unique(subtract_gr(cov_file, region_gr,keep_metadata=F))
            tiled_gr <- unlist(tile(subtract_gr(cov_file,blacklist,keep_metadata=F),width=1))

##            print(region_gr)
  ##          print(blacklist)
    ##        print(cov_file)
      ##      print(tiled_gr)
        ##    print(sum(width(subtract_gr(cov_file, blacklist, keep_metadata=F))))
            overlaps <- findOverlaps(tiled_gr,cov_file)
          ##  print(overlaps)
            tiled_gr$V4 <- cov_file$V4[overlaps@to]
            
##            print(summary(cov_file$V4))
##            print(summary(tiled_gr$V4))

        } else { cov_file <- cov_file}
        if(any(tiled_gr$V4!=0)){
            cov_vec <- tiled_gr$V4} else{ cov_vec <- rep(0,10)}
  ##      str(cov_vec)

        out <- data.frame(bam,median(cov_vec), mean(cov_vec), max(cov_vec),sum(cov_vec))
        colnames(out) <- c("Bam","Median","Mean","Max","Total")
        out2 <- as.data.frame(do.call("cbind",lapply(2:4, function(x) out[,x]/out$Total)))
##        str(out2)
        colnames(out2) <- paste0(colnames(out[2:4]),"_Norm")
        out <- cbind(out,out2)

##        str(out2)
##        str(out)
        system(paste0("rm ",out_bam))
        system(paste0("rm ",out_cov))

        return(out)}
    else{ print("Not going to estimate coverage of entire bam file here.")}
    }
            
drop_infinite <- function(df, replacement=0){

    out <- do.call(data.frame, lapply(df, function(x) replace(x, is.infinite(x), replacement)))
    rownames(out) <- rownames(df)
    return(out)}

                   

Access_pool_analysis <- function(directory,ID) {
    require(ggplot2)
    require(reshape2)
    
    
    extensions <- c("collapsed","duplex","noise","pileup","uncollapsed")
    file_list <- list.files(dir,paste0(ID,"*"))
    search_vals <- paste0(ID,"_",extensions,"*")
    file_list2 <- unique(unlist(lapply(search_vals, function(x) file_list[grep(x,file_list)])))
    file_list3 <- file_list[grep(paste0(ID,"\\."),file_list)]
    file_list <- c(file_list2,file_list3)



##    noise_list <- lapply(grep("noise",file_list,value=T), function(x) read.table(paste0(dir,x), header=T, sep='\t'))
##    names(noise_list) <- gsub(paste0(ID,"_"),"", grep("noise",file_list,value=T))
    insert_size_list <- lapply(grep("insert_size",file_list,value=T), function(x) read.table(paste0(dir,x), header=T, sep='\t',skip=10))
    names(insert_size_list) <- gsub(paste0(ID,"_"),"",grep("insert_size",file_list,value=T))
    family_size_list <- lapply(grep("family_size",file_list,value=T), function(x) read.table(paste0(dir,x), header=T, sep='\t'))
    names(family_size_list) <- gsub(paste0(ID,"_"),"",grep("family_size",file_list,value=T))
    per_target_cov_list <- lapply(grep("per_target_",file_list,value=T), function(x) read.table(paste0(dir,x), header=T, sep='\t'))
    names(per_target_cov_list) <- gsub(paste0(ID,"_"),"",grep("per_target_",file_list,value=T))
    umi_counts_list <- lapply(grep("umi_counts",file_list,value=T), function(x) read.table(paste0(dir,x), header=T, sep='\t'))
    names(umi_counts_list) <- gsub(paste0(ID,"_"),"",grep("umi_counts",file_list,value=T))



##    str(names(noise_list))
##    str(names(insert_size_list))
##    str(family_size_list)
##    str(names(per_target_cov_list))
    

    insert_size_list <- do.call("rbind",lapply(names(insert_size_list), function(x) cbind(insert_size_list[[x]],gsub("_insert_size_metrics.txt","",x))))
    colnames(insert_size_list) <- c("insert_size","reads","description")
    insert_size_list$Sample <- ID


    family_size_grouped <- family_size_list[["collapsed_grouped.family_sizes.txt"]]
    family_size_grouped <- reshape2::melt(family_size_grouped,id.vars="family_size")
    family_size_grouped$Sample <- ID
    family_size_duplex <- family_size_list[["collapsed_grouped.duplex_family_sizes.txt"]]
    family_size_duplex$Sample <- ID


    per_target_cov_list <- do.call("rbind",lapply(names(per_target_cov_list), function(x) cbind(per_target_cov_list[[x]],gsub("_per_target_coverage.txt","",x))))
    colnames(per_target_cov_list)[ncol(per_target_cov_list)] <- "description"
    per_target_cov_list$Sample <- ID
    

    umi_counts_list <- umi_counts_list[[1]]
    umi_counts_list$Sample <- ID
    

    out_list <- list(insert_size_list, family_size_grouped, family_size_duplex, per_target_cov_list, umi_counts_list)
    names(out_list) <- c("Insert_Size","Family_Size_Grouped","Family_Size_Duplex","Target_Cov","UMI_counts")
    
##    insert_fig <- ggplot(insert_size_list,aes(insert_size,reads, fill=description,alpha=0.2))
    return(out_list)

    
    }
                               

Access_pool_comparison <- function(access_pool_analysis){

    require(dplyr)
    access_pool_analysis2 <- lapply(names(access_pool_analysis[[1]]), function(x) do.call("rbind", lapply(names(access_pool_analysis), function(y) access_pool_analysis[[y]][[x]]))); names(access_pool_analysis2) <- names(access_pool_analysis[[1]])

    access_pool_analysis2$Target_Cov$Region <- paste0(access_pool_analysis2$Target_Cov$chrom,":",access_pool_analysis2$Target_Cov$start,"-",access_pool_analysis2$Target_Cov$end)
##    str(access_pool_analysis2)
    access_pool_analysis2$Target_Cov <- access_pool_analysis2$Target_Cov[c(ncol(access_pool_analysis2$Target_Cov),1:(ncol(access_pool_analysis2$Target_Cov)-1))]
##    str(access_pool_analysis2)

    access_combinations <- expand.grid(names(access_pool_analysis), names(access_pool_analysis))

    out <- list()
    for(j in names(access_pool_analysis2)){
        int <- data.frame(stringsAsFactors=F)

        for(i in 1:nrow(access_combinations)){

            sample1 <- access_combinations[i,1]
            sample2 <- access_combinations[i,2]
            sub1 <- dplyr::filter(access_pool_analysis2[[j]], Sample==sample1)
            sub2 <- dplyr::filter(access_pool_analysis2[[j]], Sample==sample2)

            merged <- merge(sub1, sub2, by=colnames(sub1)[1],all=TRUE)
            int <- rbind(int, merged)

            print(paste0("Finished ", i," of ", nrow(access_combinations)))

        }
        out[[j]] <- int
        print(paste0("Finished ",j))

    }
    return(out)}


Access_pool_comparison_plots <- function(comparison){

    require(ggplot2)
    

    }    

parse_OxBS <- function(OxBS,outfile, genome="hg19") {

    
    system(paste0("awk -v OFS='\t' '{print $1,$2, $2+2, $3}' ", OxBS,"|tail -n +2  > OxBS_int_bedgraph"))

    if(genome == "hg19"){
        system(paste0("./Andre_Software/bedGraphToBigWig OxBS_int_bedgraph /home/forbesa1/hg19_seqlengths.txt ", outfile))} else { print("Genome not defined")}

    

    print("Finished")}

group_df <- function(df, metadata, meta_col,sample_col="Sample") {
    classes <- unique(metadata[,meta_col])
    

    metadata <- metadata[which(metadata[,sample_col] %in% colnames(df)),]
    out <- as.data.frame(do.call("cbind", lapply(classes, function(x) rowMeans(df[,metadata[which(metadata[,meta_col] ==x),sample_col]]))))
##    str(out)
    colnames(out) <- classes
    return(out)}



    
rowMeans2 <- function(df){

    
    if(class(df) %in% c("data.frame","matrix")){ out <- rowMeans(df)} else{ out <- df}
##    str(out)

    return(out)}


cosine_cor <- function(df,sample_margin=c("row","column")){

    require(lsa)

    if(sample_margin=="column"){
        mat <- as.data.frame(do.call("cbind",lapply(colnames(df), function(x) unlist(lapply(colnames(df), function(y) cosine(all_peaks_all_patients_final_binary[,x], all_peaks_all_patients_final_binary[,y]))))))
        colnames(mat) <- rownames(mat) <- colnames(df)} else if(sample_margin=="row"){

mat <- as.data.frame(do.call("cbind",lapply(rownames(df), function(x) unlist(lapply(rownames(df), function(y) cosine(all_peaks_all_patients_final_binary[x,], all_peaks_all_patients_final_binary[y,]))))))
        colnames(mat) <- rownames(mat) <- rownames(df)

                                                      }


    return(mat)}
        
calc_GC <- function(gr,output=c("all","mean","median"),genome=c("hg19","hg38"),format=c("count","prob")){
    require(GenomicRanges)
    if(genome=="hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19)} else if(genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38)}
    source("~/Andre_F_functions.R")
    require(Biostrings)

    gr <- trim_invalid_intervals(gr,genome)
    sequences <- getSeq(Hsapiens,gr)

##    print(sequences[1:5])
    if(format=="count"){ freq <- Biostrings::letterFrequency(sequences,letters="GC",as.prob=FALSE)} else if(format=="prob"){ freq <- Biostrings::letterFrequency(sequences,letters="GC", as.prob=TRUE)}

    if(output=="all"){ out <- freq}  else if(output=="mean"){ out <- mean(freq)} else if(output=="median"){ out <- median(freq)}
    return(out)
    }

bigwig_callpeaks <- function(bigwig, outdir,root_name,util_path=NULL){
    source("~/Andre_F_functions.R")
    if(is.null(util_path)){
    util_path <- "/home/forbesa1/Andre_Software/bigWigToBedGraph"} else{ util_path <- util_path}
    system(paste0("mkdir ", outdir))

    bedgraph_out <- gsub("bigwig|bw","bedgraph",bigwig)
    str(bedgraph_out)


    run_cmd <- paste0(util_path," ", bigwig," ", bedgraph_out)

    print(run_cmd)
    system(run_cmd)
    system(paste0("macs2 bdgpeakcall -i ",bedgraph_out," --outdir ",outdir," --o-prefix ", root_name))
    system(paste0("rm ", bedgraph_out))
    }

smooth_gr <- function(gr, binsize=100,genome=c("hg19","hg38"),y.field="score",pad=1e3){
    require(GenomicRanges)
    require(gUtils)

        if(genome=="hg19"){
            require(BSgenome.Hsapiens.UCSC.hg19)} else if(genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38)}

        source("~/Andre_F_functions.R")

    ##        bins <- collapse_gtrack_list(tile(gr+pad,width=binsize))
    bins <- grl.unlist(tile(gr+pad,width=binsize))


        gr_out <- gr

    print(length(gr_out))
##    print(sum(width(gr_out)))

    print(length(bins))
##    print(sum(width(bins)))
    bins$score <- unlist(lapply(1:length(bins), function(x) mean(mcols(intersect_with_metadata(gr,bins[x]))[,y.field])))
    bins$score[is.na(bins$score)] <- 0

        return(bins)}

rescale_gr <- function(gr,metadata_column="score",range=c(0,1)){
    require(GenomicRanges)
    require(scales)

    mcols(gr)[,metadata_column] <- rescale(mcols(gr)[,metadata_column],to=range)

    return(gr)}

parse_Bismark_cov <- function(cov_file,outfile, genome="hg19",meta_col=4) {

    if(!grepl(".gz",cov_file)){
        system(paste0("awk -v OFS='\t' '{print $1,$2, $2+1, $",meta_col,"}' ", cov_file," > int_bedgraph"))} else{
                                                                                                               system(paste0("zcat ", cov_file," > int.cov"))
                                                                                                               system(paste0("awk -v OFS='\t' '{print $1,$2, $2+1, $",meta_col,"}' int.cov > int_bedgraph"))
                                                                                                               system("rm int.cov")}

    system("sortBed -i int_bedgraph > sorted_int_bedgraph")
    system("mv sorted_int_bedgraph int_bedgraph")


    if(genome == "hg19"){
        system(paste0("./Andre_Software/bedGraphToBigWig int_bedgraph /home/forbesa1/BergerLab_Work/Genome_files/hg19/hg19.chrom.sizes ", outfile));system("rm int_bedgraph")} else if (genome =="hg38"){ system(paste0("./Andre_Software/bedGraphToBigWig int_bedgraph /home/forbesa1/BergerLab_Work/Genome_files/hg38/hg38.chrom.sizes ", outfile)) } else{ print("Genome not defined/available")}




    print("Finished")}


generate_cox_ph_expression_basic <- function(dir,preload=NULL, input_tables=c("patient_summary","progression_events","g_rna_gene_expression","g_molecular_metadata","deceased_index"),index_table="g_molecular_metadata",index_column="biopsy_collection_date_days_from_index",general_index_root="_date_days_from_index",patient_identifier="patient_id",covariates=NULL,patient_subset=NULL,expression_df="g_rna_gene_expression",tissue_blacklist="Blood",target_gene="PLK1",gene_value_column="log2_gene_tpm_corrected",outcome=c("pfs","os")){

    require(dplyr)
    require(tempusr)
    source("Andre_F_functions.R")

    if(is.null(preload)){
        input_td <- load_tempus_data(dir, list_files=as.list(input_tables))} else if(all(input_tables %in% names(preload))){ input_td <- preload} else{ print("Missing data tables in preloaded data specified")}


    ## str(input_td)

    if(!is.null(patient_subset)){ input_td[[index_table]] <- dplyr::filter(input_td[[index_table]], get(patient_identifier) %in% patient_subset, isolate_analyte=="rna", tissue_site_canonical_name!=tissue_blacklist)} else{ dplyr::filter(input_td[[index_table]], isolate_analyte=="rna", tissue_site_canonical_name!=tissue_blacklist)}

    td_outcomes <- prepare_outcomes(input_td[[index_table]],index_column,input_td[["patient_summary"]],"last_known_followup_date_days_from_index",input_td[["deceased_index"]],paste0("deceased",general_index_root), input_td[["progression_events"]], paste0("event",general_index_root))
    ##    str(td_outcomes)
    if(!any(grepl("Gene", colnames(input_td[[expression_df]])))){
        print("Translating ENSMBL to HGNC")
        input_td[[expression_df]]$Gene <- translate_ENSMBL_to_HGNC(input_td[[expression_df]]$gene_code)
    } else { print("Not translating ENSMBL to HGNC")}
 ##  stop()
    td_outcomes <- merge(td_outcomes, dplyr::filter(input_td[[expression_df]], Gene ==target_gene),by=patient_identifier)


    tertiles <- quantile(td_outcomes[,gene_value_column],c(0.333,0.667))

    td_outcomes[,paste0(target_gene,"_Flag")] <- ifelse(td_outcomes[,gene_value_column] >= tertiles[2],paste0(target_gene,"_High"), ifelse(td_outcomes[,gene_value_column] >= tertiles[1],paste0(target_gene,"_Mid"),paste0(target_gene,"_Low")))

##    str(td_outcomes)

##    saveRDS(td_outcomes,"intermediate_pfs.rds")

## str(td_outcomes)
    if(is.null(covariates)){


    final_outcomes <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),paste0(target_gene,"_Flag"),name_outcome=toupper(outcome)))} else { stopifnot(all(covariates %in% colnames(td_outcomes[["patient_summary"]]))); final_outcomes <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),paste0(target_gene,"_Flag"),list_group_columns=list(covariates),name_outcome=toupper(outcome)))}
    ## stop()
    ##  str(final_outcomes)
    significance_continuous <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),"log2_gene_tpm_corrected",name_outcome=toupper(outcome),return_plot=FALSE)$hr)

    out <- list(td_outcomes, final_outcomes,significance_continuous)
    names(out) <- c("Intermediate_Outcomes","Final_Outcomes","Significance")
    return(out)


    }


generate_cox_ph_expression_covars <- function(dir,preload=NULL, input_tables=c("patient_summary","progression_events","g_rna_gene_expression","g_molecular_metadata","deceased_index"),index_table="g_molecular_metadata",index_column="biopsy_collection_date_days_from_index",general_index_root="_date_days_from_index",patient_identifier="patient_id",covariates=NULL,patient_subset=NULL,expression_df="g_rna_gene_expression",tissue_blacklist="Blood",target_gene="PLK1",gene_value_column="log2_gene_tpm_corrected",outcome=c("pfs","os")){

    require(dplyr)
    require(tempusr)
    source("Andre_F_functions.R")

    if(is.null(preload)){
        input_td <- load_tempus_data(dir, list_files=as.list(input_tables))} else if(all(input_tables %in% names(preload))){ input_td <- preload} else{ print("Missing data tables in preloaded data specified")}


    ## str(input_td)

    if(!is.null(patient_subset)){ input_td[[index_table]] <- dplyr::filter(input_td[[index_table]], get(patient_identifier) %in% patient_subset, isolate_analyte=="rna", tissue_site_canonical_name!=tissue_blacklist)} else{ dplyr::filter(input_td[[index_table]], isolate_analyte=="rna", tissue_site_canonical_name!=tissue_blacklist)}

    td_outcomes <- prepare_outcomes(input_td[[index_table]],index_column,input_td[["patient_summary"]],"last_known_followup_date_days_from_index",input_td[["deceased_index"]],paste0("deceased",general_index_root), input_td[["progression_events"]], paste0("event",general_index_root))
    ##    str(td_outcomes)
    if(!any(grepl("Gene", colnames(input_td[[expression_df]])))){
        print("Translating ENSMBL to HGNC")
        input_td[[expression_df]]$Gene <- translate_ENSMBL_to_HGNC(input_td[[expression_df]]$gene_code)
    } else { print("Not translating ENSMBL to HGNC")}
 ##  stop()
    td_outcomes <- merge(td_outcomes, dplyr::filter(input_td[[expression_df]], Gene ==target_gene),by=patient_identifier)


    tertiles <- quantile(td_outcomes[,gene_value_column],c(0.333,0.667))

    td_outcomes[,paste0(target_gene,"_Flag")] <- ifelse(td_outcomes[,gene_value_column] >= tertiles[2],paste0(target_gene,"_High"), ifelse(td_outcomes[,gene_value_column] >= tertiles[1],paste0(target_gene,"_Mid"),paste0(target_gene,"_Low")))

##    str(td_outcomes)

##    saveRDS(td_outcomes,"intermediate_pfs.rds")

## str(td_outcomes)
    if(is.null(covariates)){


        final_outcomes <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),paste0(target_gene,"_Flag"),name_outcome=toupper(outcome)))
        significance_continuous <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),"log2_gene_tpm_corrected",name_outcome=toupper(outcome),return_plot=FALSE)$hr)
    } else {
        covars <- unlist(strsplit(covariates,"\\+"))

        stopifnot(all(covars %in% colnames(input_td[[index_table]])))
##        print(table(covars %in% colnames(input_td[[index_table]])))
        td_merge <- merge(td_outcomes, input_td[[index_table]][c(patient_identifier,covars)],by=patient_identifier)
        saveRDS(td_merge,"test_merge.rds")

        covar_formula <- paste0(paste0(target_gene,"_Flag"),"+",paste0(covars,collapse="+"))
       ##print(covar_formula)
                                                                                                                                                                                    final_outcomes <- suppressWarnings(calc_outcomes(td_merge,paste0(outcome,"_time"),paste0(outcome,"_flag"),list_group_columns=covar_formula,name_outcome=toupper(outcome)))

                                                                                                                                                                                        significance_continuous <- suppressWarnings(calc_outcomes(td_merge,paste0(outcome,"_time"),paste0(outcome,"_flag"),paste0("log2_gene_tpm_corrected","+",paste0(covars,collapse="+")),name_outcome=toupper(outcome),return_plot=FALSE)$hr)

    }

    ## stop()
    ##  str(final_outcomes)
##    significance_continuous <- suppressWarnings(calc_outcomes(td_outcomes,paste0(outcome,"_time"),paste0(outcome,"_flag"),"log2_gene_tpm_corrected",name_outcome=toupper(outcome),return_plot=FALSE)$hr)

    out <- list(td_outcomes, final_outcomes,significance_continuous)
    names(out) <- c("Intermediate_Outcomes","Final_Outcomes","Significance")
    return(out)


    }






new_pathos_checkout <- function(cohort){
    current_checkout <- system("pathostk checkout show", intern=TRUE)[1]
##    pathostk checkout new <date> [--copy=<current_checkout_id>]


    }

new_pathos_RMD <- function(report_title,rootdir="~/",default_rmd="/Users/forbesa/ANF_default.Rmd",out_report=NULL){
    if(!is.null(out_report)){
        system(paste0("cp ", default_rmd," ",out_report))
        system(paste0("sed -i '' 's/default_title/",report_title,"/g' ",out_report))
        system(paste0("sed -i '' 's#default_rootdir#",rootdir,"#g' ",out_report))
    }

    else{ print("No output report path specified, stopping!"); stop()}}


new_prepare_hybrid_network <- function(net_name_variable, wgcna_net_rds,wgcna_beta_tom,mediation_file,BN_digraph,default_prepare_r="/Users/forbesa/p0069_prepare_network.R",out_prepare_r=NULL){
        if(!is.null(out_prepare_r)){
        system(paste0("cp ", default_prepare_r," ",out_prepare_r))
        system(paste0("sed -i '' 's/net_name_variable/",net_name_variable,"/g' ",out_prepare_r))
        system(paste0("sed -i '' 's#wgcna_net_rds#",wgcna_net_rds,"#g' ",out_prepare_r))
        system(paste0("sed -i '' 's#wgcna_beta_tom#",wgcna_beta_tom,"#g' ",out_prepare_r))
        system(paste0("sed -i '' 's#mediation_file#",mediation_file,"#g' ",out_prepare_r))
        system(paste0("sed -i '' 's#BN_digraph#",BN_digraph,"#g' ",out_prepare_r))
    }

    else{ print("No output path specified, stopping!"); stop()}}

new_create_network <- function(net_name_variable, wgcna_net_rds,wgcna_beta_tom,cm_file,BN_digraph,include_cm=TRUE,default_create_r="/Users/forbesa/p0069_create_network.R",out_create_r=NULL){

        if(!is.null(out_create_r)){
        system(paste0("cp ", default_create_r," ",out_create_r))
        system(paste0("sed -i '' 's/net_name_variable/",net_name_variable,"/g' ",out_create_r))
        system(paste0("sed -i '' 's#wgcna_net_rds#",wgcna_net_rds,"#g' ",out_create_r))
        system(paste0("sed -i '' 's#wgcna_beta_tom#",wgcna_beta_tom,"#g' ",out_create_r))
        system(paste0("sed -i '' 's#cm_file#",cm_file,"#g' ",out_create_r))
        system(paste0("sed -i '' 's#BN_digraph#",mediation_file,"#g' ",out_create_r))
    } else{ print("No output path specified, stopping!"); stop()}

    if(include_cm==TRUE){
        cm_out <- "c(TRUE,TRUE,TRUE,TRUE)"
        system(paste0("sed -i '' 's#include_cm#",cm_out,"#g' ",out_create_r))} else{
                                                                                 cm_out <- "c(TRUE,FALSE,FALSE,FALSE)"
                                                                                 system(paste0("sed -i '' 's#include_cm#",cm_out,"#g' ",out_create_r))} }


freq_table_to_matrix <- function(df, row_name, column_name,value_name="Freq",default_fill=0){
    df_mat <- matrix(default_fill, nrow=length(unique(df[,row_name])),ncol=length(unique(df[,column_name])))
    rownames(df_mat) <- unique(df[,row_name])
    colnames(df_mat) <- unique(df[,column_name])
    index <- which(df[,value_name] !=default_fill)
##    str(index)
    for(n in index){
        i=as.character(df[n,row_name])
        j=df[n,column_name]
        val=df[n,value_name]
        df_mat[i,j] <- val
    }
    return(df_mat)

}



plot_heatmap_from_freq_table <- function(df,x_column="Cancer",y_column="Disease",title=NULL,cluster_method="HC",cluster_index=c("row","column","both"),fontsize=2,freq_cutoff=0.005,grouping_column=NULL,fill_values=FALSE){

    require(dplyr)
    require(ggplot2)
    require(ComplexHeatmap)
    require(seriation)
    source("~/Andre_F_functions.R")


    if(!is.null(freq_cutoff)){
        df <- dplyr::filter(df, Freq >=max(Freq)*freq_cutoff)}
    if(is.null(grouping_column)){
    df <- df %>% group_by(get(y_column)) %>% mutate(Freq_Norm=Freq/sum(Freq)) %>% data.frame} else if (grouping_column %in% colnames(df)) {df <- df %>% group_by(get(grouping_column)) %>% mutate(Freq_Norm=Freq/sum(Freq)) %>% data.frame} else{ print("Grouping column not in provided dataframe: Stopping!"); stop()}
    ##    missing_preds <- setdiff(metadata_df[,y_column],unique(df[,y_column]))


    df[,1:2] <- apply(df[,1:2],2,as.character)
    df_mat <- freq_table_to_matrix(df,x_column,y_column,"Freq_Norm")
##    str(df_mat)
##    print(df_mat[1:5,1:5])

##    str(df)
    summary_table <- df %>% group_by(get(y_column)) %>% mutate(num_patients=sum(Freq)) %>% data.frame
    summary_table <- unique(summary_table[c(y_column,"num_patients")])
#    str(summary_table)
    col_anno <- HeatmapAnnotation(count=anno_barplot(summary_table[,"num_patients"]))

    ## str(df)
    summary_table2 <- df %>% group_by(get(x_column)) %>% mutate(num_patients2=sum(Freq)) %>% data.frame
    summary_table2 <- unique(summary_table2[c(x_column,"num_patients2")])
 ##   str(summary_table2)

    row_anno <- rowAnnotation(count=anno_barplot(summary_table2[,"num_patients2"]))

    title <- ifelse(is.null(title),"Frequency Table Heatmap",title)

    if(cluster_index =="row"){
        mat_index <- seriate(dist(df_mat), method=cluster_method)
        if(fill_values ==FALSE){
            p <- Heatmap(df_mat, name="Fraction", column_title=title, row_order=get_order(mat_index,1),cluster_columns = F,row_names_gp = gpar(font = 2,fontsize=fontsize),column_names_gp = gpar(font = 2,fontsize=12),top_annotation = col_anno, left_annotation=row_anno,rect_gp=gpar(col="white", lwd=0.5))} else if(fill_values==TRUE){
                                                  print("placeholder")}} else if(cluster_index=="column"){ if (fill_values==FALSE){ mat_index <- seriate(dist(t(df_mat)), method=cluster_method)
                                                                                                                                                                                                                                                                                                                                                                                      p <- Heatmap(df_mat, name="Fraction", column_title=title,cluster_rows = F, column_order=get_order(mat_index,1),row_names_gp = gpar(font = 2,fontsize=fontsize),column_names_gp = gpar(font = 2,fontsize=12),top_annotation = col_anno, left_annotation=row_anno,rect_gp=gpar(col="white", lwd=0.5)) } else if(fill_values==TRUE){ print("placeholder")  }} else if(cluster_index=="both"){
                                                                     mat_index1 <- seriate(dist(df_mat), method=cluster_method)
                                                                     mat_index2 <- seriate(dist(t(df_mat)), method=cluster_method)
                                                                     if(fill_values==FALSE){
                                                                                                                                                                                                              p <- Heatmap(df_mat, name="Fraction", column_title=title, row_order=get_order(mat_index1,1), column_order=get_order(mat_index2,1),row_names_gp = gpar(font = 2,fontsize=fontsize),column_names_gp = gpar(font = 2,fontsize=12),top_annotation = col_anno, left_annotation=row_anno,rect_gp=gpar(col="white", lwd=0.5))

} else if( fill_values==TRUE){print("placeholder") }}
    return(p)}

reset_ess_fonts <- function(){
    invisible(addTaskCallback(function(...) {
    if (interactive()) {
        # Remember to install crayon
        try(cat(crayon::reset("")), silent = TRUE)
    }
    TRUE
}, name = "ansi_reset"))}

get_earliest_treatment <- function(cohort,cancer_table="cancer", cancer_filter="primary",cancer_date_col="onset_date_time_days_from_index", medications_table="medications_rollup",medication_date_col="effective_date_start_days_from_index", tempus_data=NULL,biopsy_gap=90){
    require(tempusr)
    require(dplyr)
    require(reshape2)



    if(is.null(tempus_data)){
        tempus_data <- load_tempus_data(cohort, list_files=c(cancer_table,medications_table))
    } else if(c(cancer_table, medications_table) %in% names(tempus_data)){ print(paste0("Working with preloaded data for ",cohort)) }

pdx_index <- tempus_data$cancer %>% dplyr::filter(clinical_status == "primary") %>% dplyr::select(patient_id, condition_id,
           pdx_index = !!sym(cancer_date_col)) %>% unique() %>%
    # some patients have more than one pdx index
    group_by(patient_id, condition_id) %>%
    arrange(pdx_index) %>%
    dplyr::slice(1) %>%
    ungroup()

  first_medication <- tempus_data$medications_rollup |> 
    group_by(patient_id) |> 
    arrange(effective_date_start_year_indexed, !!sym(medication_date_col)) |> 
    slice(1) |> 
    dplyr::select(patient_id, 
                  med_start_index = effective_date_start_days_from_index,
                  med_start_date_index = effective_date_start_indexed,
           med_year_index = effective_date_start_year_indexed,
           drug_name=drug_class_name,
           drug_class=drug_class_group_name
           ) %>% ungroup()


    out <- list()
    out$first_medication <- dplyr::filter(first_medication, patient_id %in%  pdx_index$patient_id)
    out$patient_diagnosis <- pdx_index
    return(out)}


get_treatment_naive <- function(earliest_treatment_output){
  treatment_naive <- all_metadata |> 
    dplyr::select(patient_id, sample_id, biopsy_index = date_biopsy, biopsy_year = year_biopsy) |> 
    left_join(pdx_index, by = "patient_id") |> 
    left_join(first_medication, by = "patient_id") |> 
    mutate(treatment_naive = case_when(
      med_start_index > biopsy_index ~ TRUE, 
      (abs(biopsy_index - pdx_index) < 60) ~ TRUE, 
      TRUE ~ FALSE
    ))

  all_metadata <- all_metadata |> 
    left_join(treatment_naive |> 
                dplyr::select(sample_id, treatment_naive), 
              by = "sample_id") |> 
    mutate(treatment_naive_primary = case_when(
      treatment_naive == TRUE & primary_tumor == TRUE ~ TRUE, 
      TRUE ~ FALSE
    ))

  return(all_metadata)}

tempus_rna_expression_to_mat <- function(rna_df,gene_identifier="Gene", gene_val_column="log2_gene_tpm_corrected", patient_identifier="patient_id"){
    require(reshape2)
    require(dplyr)
    gene_check <- any(grepl("^Gene$|^gene$", colnames(rna_df)))
    if(gene_check ==TRUE){ print("Not translating Ensembl genes to Symbols")} else{
                                                                                print("Translating gene codes to symbols")
                                                                                rna_df$Gene <- translate_ENSMBL_to_HGNC(rna_df$gene_code)}

    rna_mat <- reshape2::dcast(rna_df, as.formula(paste0(gene_identifier,"~", patient_identifier)),value.var=gene_val_column, fun.aggregate=mean)
    if(length(unique(rna_mat[,gene_identifier]))==nrow(rna_mat)){
        rownames(rna_mat) <- rna_mat[,gene_identifier]
        rna_mat[,gene_identifier] <- NULL
    } else{ print("Cannot set rownames to Gene symbols")}

    return(rna_mat)}

parse_SNV_mutations <- function(df,tempus_data,table_name="g_molecular_master_file_filtered",genes_of_interest="TP53",mutation_type_blacklist=c("B","LB"),patient_subset=NULL,patient_identifier="patient_id"){

    require(dplyr)

    if(!is.null(mutation_type_blacklist)){
    sub_mut <- dplyr::filter(tempus_data[[table_name]],gene_canonical_name %in% genes_of_interest,!functional_impact %in% mutation_type_blacklist, !is.na(functional_impact))} else{ sub_mut <- dplyr::filter(tempus_data[[table_name]],gene_canonical_name %in% genes_of_interest)}

    if(!is.null(patient_subset)){ sub_mut <- dplyr::filter(sub_mut,patient_id %in% patient_subset)}

mut_mat <- as.data.frame(t(data.frame(rep("0", length(genes_of_interest)))))
rownames(mut_mat) <- "Empty"; colnames(mut_mat) <- genes_of_interest

for(i in unique(sub_mut$patient_id)){

    int <- as.data.frame(t(data.frame(rep("0", length(genes_of_interest)))))
    colnames(int) <- genes_of_interest
    sub <- dplyr::filter(sub_mut,patient_id==i)[c("result","functional_impact","gene_canonical_name")]
    int[,unique(sub$gene_canonical_name)] <- "1"

    mut_mat <- rbind(mut_mat,int)}

    rownames(mut_mat) <- c("Empty",unique(sub_mut$patient_id))
    mut_mat[is.na(mut_mat)] <- "0"
    mut_mat$patient_id <- rownames(mut_mat)


    mut_mat <- mut_mat[-1,]

    df <- dplyr::filter(df,get(patient_identifier) %in% patient_subset)
    df <- merge(df,mut_mat,by.x=patient_identifier,by.y="patient_id",all.x=T)

    }


plot_lollipop_from_tempus_mutations <- function(tempus_mmf,target_gene,output_html,title=NULL,variant_type_blacklist=NULL,assay_subset=NULL,return_df=FALSE){
    require(g3viz)
    require(dplyr)

    if(!is.null(variant_type_blacklist)){
        tempus_mmf <- dplyr::filter(tempus_mmf, functional_impact %in% variant_type_blacklist)
        if(nrow(tempus_mmf) <1){ print("Not enough mutations available after removing blacklisted types"); stop()}
}

    if(!is.null(assay_subset)){
        mutations <-  dplyr::filter(tempus_mmf, assay %in% assay_subset,
                                    gene_canonical_name==target_gene,variant_type=="Short Variant")
    } else { mutations <-  dplyr::filter(tempus_mmf,
                                         gene_canonical_name==target_gene,
                                         variant_type=="Short Variant")}
    if(nrow(mutations) <1){ print(paste0("Not enough mutations available after subsetting to ",paste0(assay_subset,collapse=","))); stop()}

    AA_translator <- data.frame(c("Ala","Arg","Asn","Asp","Cys","Gln",
                                  "Glu","Gly","His","Ile","Leu","Lys",
                                  "Met","Phe","Pro","Pyl","Ser","Sec",
                                  "Thr","Trp","Tyr","Val"),
                                c("A","R","N","D","C","Q","E","G","H",
                                  "I","L","K","M","F","P","O","S","U",
                                  "T","W","Y","V"),
                                c("Alanine","Arginine","Asparagine",
                                  "Aspartic acide","Cysteine","Glutamine",
                                  "Glumatic acid","Glycine","Histidine",
                                  "Isoleucine","Leucine","Lysine","Methionine",
                                  "Phenylalanine","Proline","Pyrolysine",
                                  "Serine","Selenocysteine","Threonine",
                                  "Tryptophan","Tyrosine","Valine"))



    colnames(AA_translator) <- c("Abbv","Letter","Full")
    rownames(AA_translator) <- AA_translator[,1]



    mutations$AA1 <- AA_translator[unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[1]))),"Letter"]

    mutations$AA2 <- AA_translator[unlist(lapply(mutations$amino_acid_change, function(x) gsub("p.","", unlist(strsplit(x,"[0-9]{1,5}"))[2]))),"Letter"]

    mutations$AA_Pos <- unlist(lapply(mutations$amino_acid_change, function(x) gsub("[A-Z]|[a-z]|\\.","",x)))

    mutations$AA_change <- paste0("p.",mutations$AA1,mutations$AA_Pos,mutations$AA2)
    mutations <- mutations[order(mutations$AA_Pos),]

    mutations$Mut_Class <- ifelse(mutations$mutation_effect=="synonymous_variant", "Silent",ifelse(mutations$mutation_effect=="missense_variant","Missense_Mutation","Other"))

    write.table(mutations,"mutations.txt", quote=F, sep='\t', row.names=F)
##    str(mutations)


    mutation_dat <- readMAF("mutations.txt", gene.symbol.col = "gene_canonical_name", variant.class.col = "Mut_Class", protein.change.col = "AA_change",sep='\t')


    plot_options <- g3Lollipop.options()
    if(!is.null(title)){ plot_options$titleText <- title}
    mutation_fig <- g3Lollipop(mutation_dat,gene.symbol=target_gene,output.filename=title,gene.symbol.col="gene_canonical_name", protein.change.col = "AA_change",btn.style="blue",plot.options=plot_options)
    htmlwidgets::saveWidget(mutation_fig,output_html)
    if(return_df){return(mutations)}
    }



call_mutations_from_bam <- function(bam, outfile, reference=NULL,regions=NULL,report_all_bases=FALSE){

    if(!file.exists(reference)){ print("Not a valid path and/or reference");stop()}
    ##This requires a working installation of bcftools and samtools
    if(!is.null(regions)){
        if(is.character(regions)){
            file_check <- file.info(regions)
            if(is.na(file_check$size)){
                coord_check <- grepl("[0-9]:[0-9]{1,9}-[0-9]{1,9}|chr[a-zA-Z0-9_]{1,100}:[0-9]{1,9}-[0-9]{1,9}",regions)
                if(coord_check){ region_flag= "-r"} else{ print("Not a valid location. Region format is \"chr:start-end");stop()}} else{ region_flag="-R"}}
        if(!report_all_bases){
        cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -vmO v -o",outfile,sep=" ")} else { cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -mO v -o",outfile,sep=" ")

                                                                                                                            }
        system(cmd)
    } else if(!report_all_bases){
        cmd <- paste("bcftools mpileup",bam,"-f",reference, "|bcftools call -vmO v -o",outfile,sep=" ")} else {cmd <- paste("bcftools mpileup",bam,region_flag,regions,"-f",reference, "|bcftools call -mO v -o",outfile,sep=" ")}
    system(cmd)
    
##    if(!grepl(".gz", outfile)){ intfile <- paste0(outfile,".gz")
##    system(paste("mv",intfile, outfile))
##      system(paste0("gunzip ",outfile))
##      } else{ outfile <- outfile}

    print(outfile)  
}


concatenate_vcfs <- function(input,sample_ID=NULL,verbose=F){
  out_vcf <- list()
  require(VariantAnnotation)
  if(is.directory(input)){
    file_list <- list.files(input,".vcf",full.names=TRUE)

    for(i in file_list){ in_vcf <- tryCatch({readVcfAsVRanges(i)}, error=function(e) { print("No variants")})
        if(!is.null(sample_ID) && sample_ID!="filename" && class(in_vcf) == "VRanges"){ sampleNames(in_vcf) <- sample_ID} else if(sample_ID=="filename" && class(in_vcf)=="VRanges"){ sampleNames(in_vcf) <- i}
        if(class(in_vcf)=="VRanges"){
    in_vcf$Sample <- unique(sampleNames(in_vcf))
    out_vcf[[as.character(in_vcf$Sample[1])]] <- in_vcf
        } else{ print("No variants")}
        if(verbose){ print(i)}}
  }
  else if( is.vector(input)){
    if(all(file.exists(input))){
      for(i in input){ in_vcf <- tryCatch({readVcfAsVRanges(i)}, error=function(e) { print("No variants")})
      if(!is.null(sample_ID) && sample_ID!="filename" && class(in_vcf) =="VRanges"){ sampleNames(in_vcf) <- sample_ID} else if(sample_ID=="filename" && class(in_vcf)=="VRanges"){ sampleNames(in_vcf) <- i}

      if(class(in_vcf)=="VRanges"){
      in_vcf$Sample <- unique(sampleNames(in_vcf))
      out_vcf[[as.character(in_vcf$Sample[1])]] <- in_vcf  } else{ print("No variants")}  
      
    if(verbose){print(i)}}}}

  out <- collapse_granges_list(out_vcf)
  return(out)
}

make_gmt_file <- function(genelist, outfile){
  if(class(genelist)=="list" && !is.null(names(genelist))){
    gmt <- unlist(lapply(names(genelist), function(x) paste0(x,"\t",paste0(genelist[[x]],collapse="\t"))))
    writeLines(gmt,outfile)
    print("Finished")
  } else{print("The provided genelist is not a *named* list object"); stop()}
}

get_gene_modules <- function(gene_info_file,gene,sep='\t',gene_column="Gene"){
  require(pathosr)
  require(tempusr)
  require(dplyr)
  ##HGNC symbol is what we're using to filter.. if you want something else, it hasn't been implemented yet
  
  if(is.character(gene_info_file)){
    print("Importing gene info file")
    gene_info <- read.table(gene_info_file,header=T, sep=sep,row.names = 1)} else{
      gene_info <- gene_info_file}
    gene_info_sub <- gene_info %>% dplyr::select(Module=module_label,Gene=hgnc_symbol,Color=module_color)
   
    if(!all(gene %in% gene_info_sub$Gene)){
      print(paste0("There is/are ",length(setdiff(gene,gene_info_sub$Gene))," gene in input missing from the gene info file: They have been removed"))
      gene <- intersect(gene,gene_info_sub$Gene)
    } else {
      gene <- gene }
   if(length(gene)<1){
     stop() 
     } else{
    module_of_interest <- unlist(lapply(gene, function(x) dplyr::filter(gene_info_sub,Gene==x)$Module))
    names(module_of_interest) <- gene
    
    return(module_of_interest)
  }
}


annotate_melanoma <- function(path, cohort_id) {
  require(tempusr)
  require(dplyr)

  # Split up by subtypes or run a model where we can include subtype as covariate
  tempus_files <- load_tempus_data(path,
                                   list_files = list("histology",
                                                     "cancer",
                                                     "g_molecular_metadata"))

  print(names(tempus_files))
  hist_rollup <- tempus_files$histology |> # one per patient
    distinct(patient_id, condition_id, hist = value_concept_canonical_name) |>
    semi_join(tempus_files$cancer |> dplyr::filter(clinical_status == "primary") |>
                dplyr::select(condition_id),
              by = "condition_id") |>
    group_by(patient_id) |>
    summarise(hist_rollup = case_when(
      any(grepl("[Aa]cral", hist)) ~ "acral",
      any(grepl("[Mm]ucosal", hist)) ~ "mucosal",
      any(grepl("[Uu]veal|orbital", hist)) ~ "uveal",
      .default = "other"
    ), .groups = "drop")

  bx_rollup <- tempus_files$g_molecular_metadata |> # one per patient
    dplyr::filter(isolate_classification == "tumor", !grepl("xF", isolate_molecular_assay)) |>
    left_join(tempusr::tissue_rollup, by = join_by(tissue_site_canonical_name)) |>
    distinct(patient_id, site = tissue_site_canonical_name, tissue_rollup) |>
    group_by(patient_id) |>
    summarize(bx_type = case_when(
      any(site %in% c(
        "Accessory sinus", "Anal canal", "Anus", "Cervix uteri", "Cheek mucosa",
        "Colon", "Esophagus", "Ethmoid sinus", "Gum", "Head, face or neck",
        "Hypopharynx", "Ileum", "Major salivary gland", "Mouth", "Nasal cavity",
        "Nasopharynx", "Overlapping lesion of rectum, anus and anal canal",
        "Palate", "Rectum", "Small intestine", "Stomach", "Submandibular gland",
        "Tongue", "Tonsil", "Transverse colon", "Vagina", "Vulva")
      ) ~ "mucosal",
      any(site %in% c("Conjunctiva", "Eye", "Eyelid", "Orbit")) ~ "uveal",
      any(tissue_rollup == "skin and external organs") ~ "cutaneous",
      all(is.na(site)) ~ "unknown",
      .default = "other"
    ), .groups = "drop")

  dx_rollup <- tempus_files$cancer |> # one per patient
    dplyr::filter(clinical_status == "primary") |>
    distinct(patient_id, rollup_organ_system_name)

  merged_subtypes <- dx_rollup |>
    left_join(hist_rollup, by = "patient_id") |>
    left_join(bx_rollup, by = "patient_id") |>
    mutate(subtype = case_when(
      hist_rollup == "acral" ~ "acral",
      hist_rollup == "mucosal" ~ "mucosal",
      rollup_organ_system_name == "Eye and orbit" ~ "uveal",
      bx_type == "uveal" ~ "uveal",
      bx_type == "mucosal" ~ "mucosal",
      rollup_organ_system_name == "Skin" ~ "cutaneous",
      bx_type == "cutaneous" ~ "cutaneous",
      .default = "other or unknown"
    )) |>
    dplyr::select(patient_id, subtype)

  return(as.data.frame(merged_subtypes))
}


quick_cohort_summary <- function(path, cohort_name,output_format=c("summary","all"),drop_versions=FALSE){
    require(tempusr)
    input_td <- load_tempus_data(path, list(c("g_molecular_metadata","metadata")))
    summary <- input_td$metadata %>% dplyr::filter(!grepl("blood|Blood", tissue_site_canonical_name)) %>% dplyr::select(patient_id,source=isolate_analyte,isolate_classification,assay=isolate_molecular_assay) %>% dplyr::filter(isolate_classification !="normal") %>% unique

    if(drop_versions==TRUE){ summary$assay <- unlist(lapply(summary$assay, function(x) unlist(strsplit(x,"\\."))[1])); summary <- unique(summary) } else{ summary <- summary}


    
    out <- summary %>% group_by(assay) %>% dplyr::summarize(Count= n()) %>% data.frame
    out$Cohort <- cohort_name

    assay_df <- unique(summary[c("source","assay")])

    rownames(assay_df) <- assay_df[,2]
    out$Source <- assay_df[out$assay,1]
    if(output_format =="summary"){
    return(out)} else if(output_format=="all") { return(list(out,summary))
    }
}

vranges_to_df_basic <- function(in_vr,ID=NULL){
    require(VariantAnnotation)
    in_vr <- unique(in_vr)
    df <- data.frame(unlist(seqnames(in_vr)),start(in_vr), end(in_vr),ref(in_vr), alt(in_vr), refDepth(in_vr),altDepth(in_vr))
    colnames(df) <- c("Chrom","Start","End","Ref","Alt","RefDepth","AltDepth")
    df$TotalDepth <- df$RefDepth+df$AltDepth
    df$AlleleFraction <- round(df$AltDepth/df$TotalDepth,3)
    if(is.null(ID)){ df$Sample <- sampleNames(in_vr)@values} else{ df$Sample <- ID}
    return(df)}

benchmark_mutation_calls <- function(in_vcf, mmf,patient_ID, assays=c("xE","RS","xT"),var_type="Short Variant",feature="TotalDepth",blacklist_regions=NULL,verbose=FALSE){
    require(dplyr)
    require(VariantAnnotation)

    source("~/Andre_F_functions_git/Andre_F_functions.R")

    if(is.character(in_vcf)){ print("Trying to import vcf");test_gr <- readVcfAsVRanges(in_vcf)} else if(class(in_vcf) %in% c("VRanges","GRanges")){ if(verbose){print("Working with the vcf as is")}; test_gr <- in_vcf} else{ print("Fatal error occurred"); stop()}


    if(!is.null(blacklist_regions)){ if(verbose){print("Removing blacklisted regions from input vcf file")}
        if(is.character(blacklist_regions)){ blacklist_regions <- bed_to_granges_dynamic(blacklist_regions)} else if(class(blacklist_regions) %in% c("GRanges","UnstitchedGPos")){
        test_gr <- setdiff_with_metadata(test_gr, blacklist_regions)} else { print("Fatal error occurred"); stop()}}


    ref_mmf <- dplyr::filter(mmf, variant_type==var_type,patient_id==patient_ID, grepl(paste0(paste0("^",assays),collapse="|"),assay))

##    str(ref_mmf)
    ref_sub <- dplyr::select(ref_mmf,chromosome,position_1,position_2,reference,alternative,somatic_germline,vaf=variant_allele_freq, coverage, gene=gene_canonical_name,nucleotide_change,AA_change=amino_acid_change, functional_impact)

    ref_sub$position_2 <- ifelse(nchar(ref_sub$position_2)==0,NA,ref_sub$position_2)
    ref_sub$position_2 <- ifelse(is.na(ref_sub$position_2), ref_sub$position_1,ref_sub$position_2)

    ref_gr <- table_to_granges(ref_sub)
    if(!is.null(blacklist_regions)){ ref_gr <- setdiff_with_metadata(ref_gr,blacklist_regions)}






    in_df <- vranges_to_df_basic(test_gr,patient_ID)

##    str(in_df)
    counter <- round(seq(1, max(in_df[,feature]),length.out=30),0)


    TP <- length(unique(ref_gr))
    FP <- 3.1e9-TP
    FN <- 0
    Precision <- TP/(TP+FP)
    Recall <-  TP/(TP+FN)

    AUPRC_out <- data.frame(TP,FP,FN,Precision,Recall,patient_ID)

    for(i in counter){
        sub_gr <- table_to_granges(dplyr::filter(in_df,get(feature)>=i))

##        print(sub_gr)
##        print(ref_gr)
        TP <- length(unique(intersect_with_metadata(unique(sub_gr), unique(ref_gr))))
        FP <- length(unique(setdiff_with_metadata(unique(sub_gr),unique(ref_gr))))
        FN <- length(setdiff_with_metadata(unique(sub_gr), gr1=unique(ref_gr)))
        Precision <- TP/(TP+FP)
        Recall <-  TP/(TP+FN)
        int <- data.frame(TP,FP,FN,Precision,Recall,patient_ID)
        AUPRC_out <- rbind(AUPRC_out, int)
       }

    return(AUPRC_out)
    }

hard_filter_vcf <- function(in_vcf, QUAL_val=10,MQBZ_val=-3,RPBZ_val=0.3,SCBZ_val=3,DP_val=NULL){
    require(dplyr)

    ##DP should be a value between 0 and 0.5 representing how much of the upper and lower ends of the distribution to remove or a value larger than 1

    in_vcf <- unique(in_vcf)
    in_vcf$index <- 1:length(in_vcf)
    if(is.null(DP_val)){
        metadata <- as.data.frame(in_vcf@elementMetadata@listData[c("QUAL","MQBZ","RPBZ","SCBZ","index")])
} else {

                                                                                                      metadata <- as.data.frame(in_vcf@elementMetadata@listData[c("QUAL","MQBZ","RPBZ","SCBZ","DP","index")])

    if(DP_val <=0.5){
                     quantiles <- quantile(metadata$DP, c(DP_val, 1-DP_val))
                     metadata <- dplyr::filter(metadata, DP >=quantiles[1] & DP <= quantiles[2])} else if (DP_val >1){ metadata <- dplyr::filter(metadata, DP > DP_val)} else{ print("Fatal error has occured with DP parameter")}

                                                                                                  }

    params <- formals()
    params <- params[2:length(params)]

##    print(params)
    metadata <- dplyr::filter(metadata, QUAL>=QUAL_val, MQBZ>=MQBZ_val,SCBZ<=SCBZ_val,abs(RPBZ)<=RPBZ_val)

    out_vcf <- in_vcf[metadata$index]

    return(out_vcf)
}




convert_genes_Mm_to_Hs <- function(gene_list, convert_from_ENSEMBL = FALSE, return_ENSEMBL = TRUE,ref_database=NULL){
    require(org.Mm.eg.db)
    require(org.Hs.eg.db)
    if(is.null(ref_database)){
        mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")} else{ mouse_human_genes <- ref_database}
##David Woods function  
#Convert input genes from ENSEMBL to Symbols

  if(convert_from_ENSEMBL == TRUE) {
    gene_list = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, 
                                      keys=gene_list, 
                                      keytype = "ENSEMBL",
                                      column = "SYMBOL")
  }
  
    output = data.frame(stringsAsFactors = F)

  for(gene in gene_list){
    class_key = mouse_human_genes  %>%  
      dplyr::filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory") %>% 
      dplyr::pull(DB.Class.Key)
    
    if(!identical(class_key, integer(0)) ){
      human_gene = mouse_human_genes %>% 
        dplyr::filter(DB.Class.Key == class_key & Common.Organism.Name== "human") %>% 
          dplyr::pull(Symbol)


if(length(human_gene) ==0){ output <- output } else{
      int <- data.frame(gene, human_gene)
      output = rbind(output, int) 
    }}}



  #Convert back to ENSEMBL
  if(return_ENSEMBL == TRUE){
    output = unname(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                          keys=output, 
                                          keytype = "SYMBOL",
                                          column = "ENSEMBL"))
  }
  return (output)
}

plot_sankey <- function(regimen_table,treatment_lines=5, num_drug_classes=6,num_drug_names=5,drop_untreated=FALSE,plot_option=c("Drug","Class","Both")){
    require(forcats)
    require(dplyr)
    require(tidyr)
    require(ggsankey)


    regimen_summary <- regimen_table %>%
        dplyr::select(patient_id,Rank=regimen_rank,Drug_Name=regimen_name,
                      Drug_Class=regimen_class,Class_Group=regimen_class_group) %>%
        group_by(patient_id) %>% mutate(Rank2=paste0(1:length(Rank),"L")) %>% ungroup %>%
        group_by(Rank2) %>% mutate(Drug_Name2 = forcats::fct_lump(Drug_Name, n = num_drug_names)) %>%
        ungroup %>% group_by(Rank2) %>%
        mutate(Class_Group2=forcats::fct_lump(Class_Group, n=num_drug_classes)) %>%
        ungroup
    regimen_summary <- regimen_summary %>% group_by(patient_id) %>%
        arrange(patient_id,Rank2) %>% ungroup %>% data.frame


    regimen_sankey <- regimen_summary %>% dplyr::select(patient_id,Treatment_Line=Rank2,Drug=Drug_Name2,Class=Class_Group2) %>% ungroup

    regimen_sankey_class <- regimen_sankey %>% tidyr::pivot_wider(id_cols=patient_id,names_from=Treatment_Line,values_from=Class)
    regimen_sankey_drug <- regimen_sankey %>% tidyr::pivot_wider(id_cols=patient_id,names_from=Treatment_Line,values_from=Drug)


    all_lines <- unique(regimen_summary$Rank2)
    treatment_cols <- intersect(paste0(1:treatment_lines,"L"),all_lines)


    regimen_sankey_class_final <- regimen_sankey_class %>% make_long(treatment_cols)
    regimen_sankey_drug_final <- regimen_sankey_drug %>% make_long(treatment_cols)



    sankey_theme <- theme(axis.text.x=element_text(face='bold',size=6,angle=15,hjust=1),axis.text.y=element_text(face='bold',size=8),strip.text = element_text(colour = "black", face = "bold",size=10),plot.title=element_text(hjust=0.5))


    if(drop_untreated){
        regimen_sankey_class_final <- regimen_sankey_class_final %>% drop_na(node)
        regimen_sankey_drug_final <- regimen_sankey_drug_final %>% drop_na(node)} else{
                                                                                    regimen_sankey_class_final <- regimen_sankey_class_final %>% replace_na(list(node="No F/U"))
                                                                                    regimen_sankey_drug_final <- regimen_sankey_drug_final %>% replace_na(list(node="No F/U"))
                                                                                    }

    if(plot_option=="Drug"){
        p2 <- ggplot(regimen_sankey_drug_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0)+sankey_theme+labs(title="Patient treatment lines by drug name")+xlab("Treatment Line")+ylab("Num. Patients")
        print(p2)} else if(plot_option=="Class"){
                     p <- ggplot(regimen_sankey_class_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0)+sankey_theme+labs(title="Patient treatment lines by drug class")+xlab("Treatment Line")+ylab("Num. Patients")
                     print(p)
} else if(plot_option=="Both"){

    p <- ggplot(regimen_sankey_class_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0.25,position=position_nudge(x=0.1))+sankey_theme+labs(title="Patient treatment lines by drug class")+xlab("Treatment Line")+ylab("Num. Patients")
                     print(p)

    p2 <- ggplot(regimen_sankey_drug_final,aes(x=x, next_x=next_x, node=node,next_node=next_node,label=node,fill=as.factor(node)))+geom_sankey(show.legend=F,node.color=1)+geom_sankey_text(size = 3, color = "black", hjust = 0.25,position=position_nudge(x=0.1))+sankey_theme+labs(title="Patient treatment lines by drug name",)+xlab("Treatment Line")+ylab("Num. Patients")
    print(p2)
}

    }



