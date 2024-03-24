##' Rex plots
##'
##'
##'
##'
##'
##'
##'
plotVolcano <- function(diff_params,
    nrow = 5,
    quantity = "TRE") {
    if (quantity == "TRE") {
        TRE <- diff_params@Rex.estimates[, grep("TRE", colnames(diff_params@Rex.estimates))]
        probs <- diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))]
        timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(TRE)))



        df_volcanoe <- data.frame(
            TRE = unlist(TRE),
            probs = unlist(probs),
            residues = rep(diff_params@Rex.probs$Residues,
                times = length(timepoints)
            ),
            timepoints = rep(timepoints,
                each = length(diff_params@Rex.probs$Residues)
            )
        )


        gg1 <- df_volcanoe %>%
            ggplot(aes(x = TRE, y = probs, col = TRE < 0)) +
            geom_point(size = 2) +
            theme_bw() +
            facet_wrap(. ~ timepoints, nrow = nrow) +
            xlim(c(-max(abs(df_volcanoe$TRE)), max(abs(df_volcanoe$TRE)))) +
            xlab("Effect Size (TRE)") +
            ylab("Probability") +
            scale_color_manual(
                values = alpha(c("darkred", "darkblue"), 0.5),
                labels = c("Deproctetion", "Protection")
            ) +
            labs(color = "") +
            guides(
                alpha = guide_legend(override.aes = list(size = 5)),
                color = guide_legend(override.aes = list(size = 5))
            ) +
            theme(
                text = element_text(size = 20),
                panel.spacing = unit(1, "lines"),
                strip.background = element_rect(fill = "steelblue")
            )
    } else if (quantity == "ARE") {
        ARE <- diff_params@Rex.estimates[, grep("ARE", colnames(diff_params@Rex.estimates))]
        probs <- diff_params@Rex.probs[, -c(1, ncol(diff_params@Rex.probs))]
        timepoints <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(ARE)))



        df_volcanoe <- data.frame(
            ARE = unlist(ARE),
            probs = unlist(probs),
            residues = rep(diff_params@Rex.probs$Residues,
                times = length(timepoints)
            ),
            timepoints = rep(timepoints,
                each = length(diff_params@Rex.probs$Residues)
            )
        )

        gg1 <- df_volcanoe %>%
            ggplot(aes(x = ARE, y = probs, col = ARE < 0)) +
            geom_point(size = 2) +
            theme_bw() +
            facet_wrap(. ~ timepoints, nrow = nrow) +
            xlim(c(-max(abs(df_volcanoe$ARE)), max(abs(df_volcanoe$ARE)))) +
            xlab("Effect Size (ARE)") +
            ylab("Probability") +
            scale_color_manual(
                values = alpha(c("darkred", "darkblue"), 0.5),
                labels = c("Deproctetion", "Protection")
            ) +
            labs(color = "") +
            guides(
                alpha = guide_legend(override.aes = list(size = 5)),
                color = guide_legend(override.aes = list(size = 5))
            ) +
            theme(
                text = element_text(size = 20),
                panel.spacing = unit(1, "lines"),
                strip.background = element_rect(fill = "steelblue")
            )
    }



    return(gg1)
}
