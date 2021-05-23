
volcano_HH <- function(in_test = x1, color_range = c(x2,x3,x4), mycolors = c(x5,x6,x7), variable_label = x8, in_width = x9, in_height = x10){
  
  list.of.packages <- c("dplyr", "readxl", "ggplot2", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

  ##################################[ Code ]##################################

  ####################[ Function ]####################

  ttest <- function(df, grp1, grp2) {
    x = df[grp1]
    y = df[grp2]
    x = as.numeric(x)
    y = as.numeric(y)
    results = t.test(x, y)
    results$p.value
  }


  utest <- function(df, grp1, grp2) {
    x = df[grp1]
    y = df[grp2]
    x = as.numeric(x)
    y = as.numeric(y)
    results = wilcox.test(x, y,
                          alternative = c("two.sided"),
                          mu = 0,
                          con.int = FALSE,
                          conf.level = 0.95)
    results$p.value
  }


  check_include <- function(tar_data, find_factor) {
    x <- tar_data
    y <- find_factor
    z <- grepl(y, names(results_spl))
    for(t in 1:length(z)) {
      if(z[t] == TRUE) {
        return(1)
        break
      }
    }
    return(0)
  }

  ####################[ check file extension ]####################

  file_extension <- unlist(strsplit(in_file, split="[.]"))
  if(file_extension[2] == "csv") {
    df <- read.csv(in_file, check.names = FALSE)
  } else {
    df <- read_excel(in_file)
  }

  ################################################################


  spl <- split(df, df[,2])


  ####################[ total for ]####################

  for(i in 1:length(spl)) {
    for(j in 1:length(spl)) {
      if(i==j) {

      } else {


        df_bind <- as.data.frame(rbind(spl[[i]],spl[[j]]))
        dimnames(df_bind)[[1]] <- df_bind[,1]
        df_bind_minus <- df_bind[,-1:-2]
        df_t <- t(df_bind_minus)



        dr1 <- nrow(spl[[i]])
        dr2 <- nrow(spl[[j]])





        ########################[ P-value ]########################

        if(in_test=="ttest") {

          ####################[ ttest ]####################
          rawpvalue = apply(df_t, 1, ttest, grp1 = c(1:dr1), grp2 = c((dr1+1):(dr1+dr2)))
          #hist(rawpvalue)

        } else {

          ####################[ utest ]####################
          rawpvalue = apply(df_t, 1, utest, grp1 = c(1:dr1), grp2 = c((dr1+1):(dr1+dr2)))
          #hist(rawpvalue)

        }

        ###########################################################












        #########################[ Fold change ]##########################

        df_t = log2(df_t)
        control = apply(df_t[,1:dr1], 1, mean)
        test = apply(df_t[, (dr1+1):(dr1+dr2)], 1, mean)

        foldchange <- control - test
        #hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

        results = cbind(foldchange, rawpvalue)
        results = as.data.frame(results)
        results$probename <- rownames(results)

        ##################################################################









        ####################[ Color Up/Down/No ]####################

        c(pval <- color_range[1], fold1 <- color_range[2], fold2 <- color_range[3])

        results$udn <- "NO"
        results$udn[results$foldchange >= fold2 & results$rawpvalue < pval] <- "UP"
        results$udn[results$foldchange <= fold1 & results$rawpvalue < pval] <- "DOWN"

        results$delabel <- NA

        results_spl <- split(results, results$udn)

        names(mycolors) <- c("DOWN", "UP", "NO")

        ############################################################












        ####################[ Top10 & Low10 / Probename ]####################



        if(check_include(names(results_spl),"DOWN") == 1) {
          if(check_include(names(results_spl),"UP") == 1) {
            low_10 <- tail(arrange(results_spl[["DOWN"]], desc(foldchange)), 10)
            top_10 <- head(arrange(results_spl[["UP"]], desc(foldchange)), 10)

            for(k in 1:nrow(results)) {
              for(l in 1:nrow(top_10)) {
                if(results$probename[k] == top_10[l,3]) {
                  results$delabel[k] <- top_10[[l,3]]
                }
              }
            }
            for(m in 1:nrow(results)) {
              for(n in 1:nrow(low_10)) {
                if(results$probename[m] == low_10[n,3]) {
                  results$delabel[m] <- low_10[[n,3]]
                }
              }
            }

          } else {
            low_10 <- tail(arrange(results_spl[["DOWN"]], desc(foldchange)), 10)

            for(m in 1:nrow(results)) {
              for(n in 1:nrow(low_10)) {
                if(results$probename[m] == low_10[n,3]) {
                  results$delabel[m] <- low_10[[n,3]]
                }
              }
            }
          }
        } else {
          if(check_include(names(results_spl),"UP") == 1) {
            top_10 <- head(arrange(results_spl[["UP"]], desc(foldchange)), 10)

            for(k in 1:nrow(results)) {
              for(l in 1:nrow(top_10)) {
                if(results$probename[k] == top_10[l,3]) {
                  results$delabel[k] <- top_10[[l,3]]
                }
              }
            }
          }
        }



        #####################################################################














        ##########################[ Volcano plot ]###########################

        volcano = ggplot(data = results,
                         aes(x = foldchange, y = -1*log10(rawpvalue),
                             col = udn, label = delabel))

        cart1 <- ifelse(
          abs(max(results$foldchange)) - abs(min(results$foldchange)) >= 0 ,
          abs(max(results$foldchange)),
          abs(min(results$foldchange)))



        volcano = volcano +
          geom_point(size = 2) +
          geom_vline(xintercept=c(fold1, fold2), col="red") +
          geom_hline(yintercept=-log10(pval), col="red") +
          coord_cartesian(xlim = c(-cart1, cart1)) +
          scale_colour_manual(values = mycolors) +
          theme_minimal() +
          labs(title="Volcano plot") +
          theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, face = "bold")) +
          theme(axis.title = element_text(size= 20)) +
          theme(axis.text = element_text(size= 14)) +
          guides(color=guide_legend(title = NULL)) +
          xlab(bquote(Log[2]~"fold change"))+
          ylab(bquote(-Log[10]~"P"))





        if (variable_label == "NO") {
          volcano
        } else {
          volcano = volcano +
            geom_label_repel(min.segment.length = 0, seed = 42, box.padding = 0.5,
                             max.time = 1, max.iter = 1e5, point.padding = 0,
                             size = 5)
        }


        #####################################################################















        ##################################[ Save ]##################################

        setwd(ex_path)
        ifelse(dir.exists("volcano_plot"), F, dir.create("volcano_plot"))


        setwd(paste0(ex_path, "/volcano_plot"))
        write.csv(results[,-3:-5],
                  file = paste0("(",spl[[i]][1,2],")_per_(",spl[[j]][1,2],")_", in_test, "_.csv"),
                  row.name = TRUE)
        png(filename= paste0("(", spl[[i]][1,2],")_per_(",spl[[j]][1,2],")_", in_test, "_.png"),
            width=in_width, height=in_height,unit="px")
        plot(volcano)
        dev.off()


        setwd(ex_path)


        ############################################################################

      }
    }
  }
}
