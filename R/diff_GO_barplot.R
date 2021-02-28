#' This function is used to plot different GO terms between two groups
#' @details barplot part which exhibit the topn different expressed genesets in two group
#' @param up_path file(.csv) of up pathways in GO analyse
#' @param down_path file(.csv) of down pathways in GO analyse
#' @param topn number of top pathways to show
#' @param width the width of the figure
#' @param height the height of the figure
#' @param text.size the size of the text
#' @param filename the name of the generated file
#'
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' \dontrun{
#' up_path="../GO_up_GO_UP_result.csv"
#' down_path="../GO_down_GO_UP_result.csv"
#' GO_diff_barplot(up_path = up_path,down_path = down_path,topn = 10,
#                filename = "barplot",text.size = 2,width = 8,height = 6)
#' }
#'
GO_diff_barplot = function(up_path=NULL,down_path=NULL,topn=10,width=8,height=6,text.size=1.5,
                         filename="barplot"){
  p.adjust <- group <- pathname <- NULL
  UP <- read.csv(up_path)
  UP$p.adjust <- -log10(UP$p.adjust)
  DOWN <- read.csv(down_path)
  DOWN$p.adjust <- -log10(DOWN$p.adjust)
  UP$group <- rep("UP",nrow(UP))
  DOWN$group <- rep("DOWN",nrow(DOWN))
  UP <- UP[1:topn,]
  DOWN <- DOWN[1:topn,]

  UP_DOWN <- rbind(UP,DOWN)
  UP_DOWN$p.adjust[which(UP_DOWN$group=="DOWN")] <- 0-UP_DOWN$p.adjust[which(UP_DOWN$group=="DOWN")]
  UP_DOWN$pathname <- UP_DOWN$Description
  UP_DOWN <- UP_DOWN[order(UP_DOWN$p.adjust,decreasing = F),]
  UP_DOWN$group <- factor(UP_DOWN$group)
  UP_DOWN$pathname <- factor(UP_DOWN$pathname,levels = UP_DOWN$pathname)

  g1 <- ggplot(UP_DOWN,aes(pathname,p.adjust,fill=group))+geom_bar(stat = "identity")+
    scale_fill_manual(values = c("lightskyblue","hotpink2"))+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
    theme(axis.title.y = element_blank())+
    theme(axis.ticks.y =  element_blank())+
    theme(axis.line.y = element_blank())+
    theme(axis.text.y = element_blank())+
    ylab("Pathname")+
    coord_flip()

  g1 <- g1+
    geom_text(aes(y=c(rep(0.2,nrow(UP_DOWN))),label=c(as.character(pathname[1:topn]),
                                                      rep("",topn)),
                  vjust = 0.5,hjust=0),size=text.size)+
    geom_text(aes(y=c(rep(-0.2,nrow(UP_DOWN))),label=c(rep("",topn),
                                                       as.character(pathname[(topn+1):nrow(UP_DOWN)])),
                  vjust = 0.5,hjust=1),size=text.size)

  pdf(paste0(filename,"_GO_path_diff.pdf"),width =width, height = height)
  print(g1)
  dev.off()
}

