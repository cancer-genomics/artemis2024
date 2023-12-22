
#' @export
get_cv_preds<-function(features,model) {
	features <- features %>% dplyr::mutate(rowIndex = 1:n())
	ids <- features %>% select(id, rowIndex,type)
	preds <- model$pred
	preds <- preds %>% dplyr::group_by(rowIndex) %>% dplyr::summarize(score = mean(cancer))
	preds <- inner_join(ids, preds, by="rowIndex")
	preds <- preds %>% select(-rowIndex)
	preds
}

#' @export
get_test_preds<-function(features,model) {
  preds<-predict(model,features,type="prob")
  res<-tibble(id=test$id,type=test$type,score=preds$cancer)
  res
}


###got this from rlucas
#' @export
roc_info <- function(obs, score) {
    roc <- pROC::roc
    roc <- suppressMessages(roc(obs, score,
                                levels=c("healthy", "cancer"),
                                ci=TRUE))
    list(sens = rev(roc[["sensitivities"]]),
         spec = rev(roc[["specificities"]]),
         auc = roc$auc,
         lower = roc$ci[1],
         upper = roc$ci[3],
         thresholds=rev(roc[["thresholds"]]))
}



#' @export
format_roc<-function(type,scores,name,l=c("healthy","cancer")) {
  roc_val<-pROC::roc(type,scores,ci=TRUE,levels=l)
  roc_type<-paste0(name,": ",round(roc_val$auc[1],2)," (",round(roc_val$ci[1],2),"-",round(roc_val$ci[3],2),")")
  res<-tibble(type=roc_type,sens=roc_val$sensitivities,spec=roc_val$specificities)
  res
}

#' @export
plot_roc<-function(results) {
  colors <- brewer.pal(length(unique(results$type)), "Dark2")

  roc_colors <- colors
  results<-results %>% dplyr::group_by(type,spec)%>% dplyr::mutate(sens=sort(sens))

  B <- results %>%
    ggplot(aes(spec, sens, group=type,color=type)) +
    geom_vline(xintercept=0.80,
               color="gray80", size=0.5, linetype="dashed") +
    geom_line(aes(color=type), size=1.1) +
    scale_x_reverse(expand=c(0, 0.01),
                    breaks=c(0, 0.25, 0.5, 0.80, 1),
                    labels=as.character(
                        c("0", ".25", ".50", ".80", "1.0"))) +
    scale_y_continuous(expand=c(0, 0.01),
                       breaks=c(0, 0.25, 0.5, 0.75, 1),
                       labels=as.character(
                           c("0", ".25", ".50", ".75", "1.0"))) +
    scale_color_manual(values=roc_colors) +
    theme_classic(base_size=20) +
    theme(panel.grid=element_blank(),
          legend.position=c(0.6, 0.2),
          aspect.ratio=0.8,
          legend.text.align=1,
          legend.title=element_text(size=16)) +
    xlab("Specificity") + ylab("Sensitivity") +
    guides(color=guide_legend(title="AUC: (95% CI)", hjust=1))
B

}

#' @export
rocstats <- function(obs, score) {
    roc <- pROC::roc
    roc <- suppressMessages(roc(obs, score,
                                levels=c("healthy", "cancer"),
                                ci=TRUE))
    list(sens = rev(roc[["sensitivities"]]),
         spec = rev(roc[["specificities"]]),
         thresh = rev(roc$threshold),
         auc = roc$ci[2],
         lower = roc$ci[1],
         upper = roc$ci[3])
}


#' @export
get_score_threshold<-function(data,spec) {
  r<-rocstats(data$type,data$score)
  m<-r$thresh[length(r$spec[round(r$spec,2)>=spec])]
  m
}

