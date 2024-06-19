plot_tsne <- function(df,title=NULL,color="Cluster",label="Cluster",viz=c("basic","medium","enhanced"),shape=NULL,alpha=0.5,size=3,font_size=3,legend_col=2){
    require(ggplot2)
#    require(ggforce)
    require(ggrepel)
    require(dplyr)

    df <- as.data.frame(df)
    if(all(class(unlist(df[,color])) %in% c("factor","integer"))){ df[,color] <- as.character(df[,color])}
    if(color != label){
        if(is.null(shape)){
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
        p <- ggplot(df, aes(Dim1,Dim2, fill=!!sym(color),label=!!sym(label)))+geom_point(size=size, alpha=alpha,shape=21,color="gray50")+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=legend_col,override.aes=list(shape=21)))+geom_mark_rect(aes(label=!!sym(label),fill=!!sym(color)),show.legend=FALSE,label.buffer=unit(1,"mm"),expand=unit(3,"mm"),label.fontsize=10,con.type="straight")} else if(viz=="basic"){
            p <- ggplot(df, aes(Dim1,Dim2, fill=!!sym(color),label=!!sym(label)))+geom_point(size=size,shape=21,color="gray50",stroke=0.1,alpha=alpha)+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=legend_col,override.aes=list(shape=21)))+geom_text(fontface="bold",size=font_size,check_overlap=TRUE,show.legend=FALSE)} else if(viz=="medium"){
                                                                                                                                                                                                                                                                                                                                                                                                                                                     if(is.null(shape)){ p <- ggplot(df, aes(Dim1,Dim2,fill=!!sym(color)))+geom_point(size=size,alpha=alpha,shape=21, color="gray50")+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(fill=guide_legend(ncol=legend_col))+geom_text_repel(data=centroid_df, aes(label=!!sym(label),color=!!sym(color)), bg.color="white", bg.r=0.05,segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE,nudge_y=0.04*(range(df$Dim2)[2]-range(df$Dim2)[1]))
                                                                                                                                                                                                                                                                                                                                                                                                                                                         }
                else if(!is.null(shape)){p <- ggplot(df, aes(Dim1,Dim2,label=!!sym(label),fill=!!sym(color)))+geom_point(aes(fill=!!sym(color),shape=!!sym(shape)),size=size,alpha=alpha,color="gray50")+theme_classic()+theme(axis.text.x=element_text(face='bold',size=15),axis.text.y=element_text(face='bold',size=15))+labs(title=title)+guides(color=FALSE,fill=guide_legend(ncol=legend_col,override.aes=list(shape=21)))+scale_shape_manual(values=c(21:25,12:14))+geom_text_repel(data=centroid_df, aes(label=!!sym(label),color=!!sym(color)), segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE,nudge_y=0.04*(range(df$Dim2)[2]-range(df$Dim2)[1]))


                    ###+geom_text_repel(data=centroid_df, aes_string(label=label), segment.size=0.2,fontface="bold",size=font_size,show.legend=FALSE)
                }
            }




    return(p)
}


TSNE_scale <- function(color_list=NULL,extension=".txt",scales=c("fill","color"),...){
    require(ggplot2)
##    print(color_list)
    if(is.null(color_list)==TRUE){
##        print("1")
        cols <- readRDS("~/ANF_color_scale.rds")} else if(is.character(color_list) ==TRUE & !is.vector(color_list)) { if(extension == ".rds") { print("reading RDS"); cols <- readRDS(color_list)} else if(extension==".txt"){ print("reading txt");cols <- read.table(color_list,sep='\t',stringsAsFactors=F,header=T,comment.char="$")}
                                                                                                      if(is.data.frame(cols) ==TRUE){ df <- cols; cols <- df[,2]; names(cols) <- df[,1]} else{ cols <- cols}} else if(is.vector(color_list)==TRUE){ cols <- color_list}


    g <- ggplot2:::manual_scale(scales, values=cols)
    return(g)
}
