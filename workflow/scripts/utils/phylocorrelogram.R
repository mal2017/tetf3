


.crlg_to_plot_df <- \(x) set_names(as_tibble(x$res),c("x","ci.upper","ci.lower","y"))


plot_crlg <- function(crlg) {
  h0 <- -1/(crlg$n - 1) # see the phylosignal paper/code for the expected null value calculation
  df <- .crlg_to_plot_df(crlg)
  ggplot(df,aes(x, y=y)) +
    geom_path() +
    geom_path(aes(y=ci.upper), linetype="dotted") +
    geom_path(aes(y=ci.lower), linetype="dotted") +
    xlab("patristic distance") + ylab("autocorrelation") +
    geom_hline(yintercept = h0, color="gray")
}
