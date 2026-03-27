#library(trace)
library(Fragman) # This one has better documentation so we adopt it right now

# We get a list of all the fsa files from Fra
DIR <- '/mnt/trcanmed/snaketree/prj/Sanger_traces/local/share/data/run_0326'
fsas <- list.files(DIR, full.names = TRUE, pattern = ".fsa")

### Code for manual inspection of fsa files to determine which is the ladder channel for each of them 
#d_fsas <- read_fsa(fsas) # this is a trace function that under the hood uses read.abif from the seqinr package,
# also read_fsa from Fragman uses it 
d0 <- seqinr::read.abif(fsas[10], verbose=T) # same as d_fsas[[1]]$fsa
d1 <- d0$Data$DATA.1
d2 <- d0$Data$DATA.2
d3 <- d0$Data$DATA.3
d4 <- d0$Data$DATA.4

n <- length(d1)
pd <- data.frame(intensity=c(d1, d2, d3, d4), fluo=c(rep('fl', n), rep('joe', n), rep('tamra', n), rep('cxr',n)),
                 x=rep(seq(1, n), 4))
ggplot(data=pd, aes(y=intensity, x=x, color=fluo))+geom_line()
# this is plotting data from the first fsa only

# Note that there are more d0$Data$DATA.1-4 (also DATA.5-8 in the first example) but only the first four with 'sensible' numbers

# These are directly the raw data stored in the files. Which fluorofore is there / where is the ladder is not clear
# it was in the fourth channel in the first data, I would check looking at the data directly where is the ladder.
# We need a ladder for each fsa file to map these signals (all channels of a single fsa file share the length
# of the signal in different channels) to basepairs
# 
# From now on fragman code
fd <- storing.inds(DIR, lets.pullup=F)
# we change the default lets.pullup and set is as false
#lets.pullup
#A FALSE/TRUE value indicating if data should be treated for noise from channel to channel known as pull up or pull down peaks since wavelengths where the dyes are read usually overlap (blue->green->yellow->red->orange. The dafault is TRUE
# since we saw some of the previous signals to be 'strangely' modified 
#https://docs.google.com/presentation/d/1PMQWcPEmx1dVWgOFLy_crpC2k-jW654QFYbkU9bpXJ0/edit?usp=sharing
# slide 5-6-7
#in any case we shouldn't have fragments/peaks of the same lengths in the same mix

# Other steps: a smoothing Fourier transform and a correction for saturated peaks (code still not analyzed)

# These are the lengths of the ladder Fra used in the first experiment (I expect it to be the same but we need to check)
ladder <-  unique(c(seq(60, 200, by=20), seq(200, 500, by=25), seq(500,600, by=50)))
# here we need to manually put channel.ladder=4 from the previous visual inspection
# trace has a trace:::find_ladders(d_fsas, config, show_progress_bar = TRUE) # questo aggiunge $trace_bp_df che è dove sono cercati i picchi.
# function but I'd start simple
ladder.info.attach(stored=fd, ladder=ladder, channel.ladder=4)

# save this plot!
# this add info to the fd object and under the hood
# does a linear regression of the ladder signal that predicts the molecular weights from the absolute x coords.
# then this trained model is used to assign mw also to other channels (function) 
# if channel.ladder is not specified this function assume it's the last channel with data (would have been right for our first files, not sure if it's always true by construction).



# In the objects I found some info on the wavelengths for different channels but they were never == to the 
# ones in Fra notes:
# $Data$DyeW.1
# [1] 540   +20 rispetto a fine flues
# 
# $Data$DyeW.2
# [1] 568 +20 rispetto a fine joe
# 
# $Data$DyeW.3
# [1] 595 +18 rispetto a fine tamra
# 
# $Data$DyeW.4
# [1] 615 +8 rispetto a fine ladder
# 

#
# Now fsa stored is a list-like object with one element for fsa file
# elements are matrixes with the corrected signals  and rows are lenghts of the runs

# overview2 plots the data mapping those lengths to m.w. using the previous models
# it allows clicking on the trace to select peaks than are then returned:
# xlim here was focused knowing where the signal was expected for different channels
overview2(my.inds=fd[1], channel = 1:3, ladder=ladder, lwd=2, channel.ladder=4, xlim=c(80, 200), my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))
# in this case we have three automatically defined peaks only for channel_1
# parameters not all investigated
# launch this instruction and click on all the peaks you see, then press esc, this creates positions that will then be scored
peaks <- locator(type='p', pch=20, col='red')$x
# locator is broken? it gets smaller x coordinates today, it was working when I worked on it previously
# resize of plot window in rstudio server? Yes, it gets easily broken, issues
# can be detected if red dots appear not where one clicked when using locator
#### mhhh.
#peaks <- c(99) # we could skip it and put the expected lengths of polyG on the genome/from Fra
# caution, it's really very sensitive (99 works, 100 doesn't)
res <- score.markers(my.inds=fd[1], channel = 1, ladder=ladder, panel=peaks, electro=T, my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))

# res has in pos the absolute x, in hei the stregth of the signal in that position and in wei the 
# predicted (using the ladder) mol weights for them.

# get.scores(res) # I do not think we need this.

## do we dig in overview2 to extract bp - signal and work outside of R? I do not know anything about peak detection algorithms.
# we want:
#new.whole.data[[h]] <- list(xx = newxx, yy = newyy)

###########
getdata <- function (my.inds, channel = 1, ladder, xlim = NULL, ylim = NULL, 
          n.inds = NULL, channel.ladder = NULL, ploidy = 2, method = "iter2", 
          init.thresh = NULL, ladd.init.thresh = 200, lwd = 0.25, warn = TRUE, 
          min.panel = 100, suggested = TRUE, env = parent.frame(), 
          my.palette = NULL, verbose = TRUE) 
{
    ### indolore se non ci sono limiti li setta dai dati ####V
    if (length(channel) > 1 & is.null(ylim)) {
    limosna <- c(min(unlist(lapply(my.inds, min, na.rm = TRUE)), 
                     na.rm = TRUE), max(unlist(lapply(my.inds, max, na.rm = TRUE)), 
                                        na.rm = TRUE))
  }
  else {
    limosna <- ylim
  }
  ############################################################V 
  suggested.list <- list()
  counter <- 0
  ############### itera su tutti i canali, ma in teoria noi gliene passiamo soltanto 1
  for (hhh in channel) {
    counter <- counter + 1
    cols <- hhh
    dev = 50
    oldw <- getOption("warn")
    options(warn = -1)
    # noi non lo cambiamo il method quindi ok
    if (method == "ci") {
      print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
    }
    ### A parte gli error my friend blocco indolore che controlla se il canale del ladder indicato è sano e se non è indicato
    #prende l'ultimo come default
    if (is.null(channel.ladder)) {
      channel.ladder <- dim(my.inds[[1]])[2]
    }
    else {
      channel.ladder <- channel.ladder
    }
    if (dim(my.inds[[1]])[2] < channel.ladder) {
      print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channel/colors"))
      stop
    }
    ############################################################
    #### se non si passa un init.thresh e suggested è true se lo prende come ultimo percentile di tutti i dati di questo canale
    if (is.null(init.thresh) & suggested) {
      listaaaa <- do.call("cbind", lapply(my.inds, function(x) {
        y <- x[, cols]
        return(y)
      }))
      init.thresh <- quantile(listaaaa, 0.99)
    }
    ##### ok se un odice di voler vedere meno fsa di quelli totali lo fa, setta gli xlim in base al ladder se non sono definiti come argomento
    #gestione progress bar
    if (is.null(n.inds)) {
      n.inds <- c(1:length(my.inds))
    }
    else {
      n.inds <- n.inds
    }
    if (is.null(xlim)) {
      xlim <- c(min(ladder), max(ladder))
    }
    else {
      xlim <- xlim
    }
    count <- 0
    tot <- length(n.inds)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
    ###############################V
    # rimappa se si usano meno fsa dei totali il tutto togliendo quelli non utilizzati
    my.inds2 <- list(NA)
    for (i in 1:length(n.inds)) {
      v1 <- n.inds[i]
      my.inds2[[i]] <- my.inds[[v1]]
      names(my.inds2)[i] <- names(my.inds)[i]
    }
    my.inds <- my.inds2
    ################ ???
    ncfp <- c("COL1", "COL2", "COL3", "COL4", "COL5")
    
    #################V
    #####################colori plot
    if (!is.null(my.palette)) {
      cfp <- rep(my.palette, 100)
    }
    else {
      cfp <- c("cornflowerblue", "chartreuse4", "gold2", 
               "red", "orange", "purple")
    }
    col.list <- list(NA)
    att1 <- numeric()
    list.data <- list(NA)    ###############################V
    # list.data.covarrubias è popolato nell'env quando si setta il ladder prima, se lo recupera in list.data
    # o lo risetta nello stesso modo in cui lo fa ladder.info.attach (se uno lo volesse skippare, credo)
    if (exists("list.data.covarrubias")) {
      list.data <- env$list.data.covarrubias
    }
    else {
      list.ladders <- lapply(my.inds, function(x) {
        y <- x[, channel.ladder]
        return(y)
      })
      list.data <- lapply(list.ladders, find.ladder, ladder = ladder, 
                          draw = F, dev = dev, warn = warn, method = method, 
                          init.thresh = ladd.init.thresh)
    }
    # è usato questo
    # predice in base a x assoluta delle posizioni dei picchi i loro pesi molecolari con sta polinomiale fino x⁵
    # questo si rifà tutte le volte XXX FIXME
    list.models <- lapply(list.data, function(da) {
      y <- da[[3]] # 2 is hei 3 is wei?
      x <- da[[1]]
      mod <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5), 
                data = da)
      return(mod)
    })
    # questo mai
    list.models.inv <- lapply(list.data, function(da) {
      x <- da[[3]]
      y <- da[[1]]
      mod <- lm(y ~ x, data = da)
      return(mod)
    })
    # i modelli girano sul canale del ladder ad ogni chiamata di overview per i diversi canali, vabbè.
    # xx per ogni oggetto fsa ha 1:sua lunghezza run, lavoriamo con 1 sola cols per volta dato il for su tutte a riga 118
    xx <- lapply(my.inds2, function(x, cols) {
      1:length(x[, cols])
    }, cols = cols)
    newxx <- numeric()
    newyy <- numeric()
    new.whole.data <- list(NA)
    # per tutte le x
    for (h in 1:length(xx)) {
      h1 <- n.inds[h]
      count <- count + 1
      # predice il MW (perchè in un loop?)
      newxx <- as.vector(try(predict(list.models[[h1]], 
                                     newdata = data.frame(x = xx[[h]])), silent = TRUE))
      newyy <- my.inds2[[h]][, cols]
      ## questo è l'oggetto che vogliamo
      new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
      setTxtProgressBar(pb, (count/tot))
    }
    # rimappa i wei in base all'xlim inferiore per plottare a partire da dove chiediamo (?)
    common <- lapply(list.data, function(x, xlim) {
      mins <- abs(x$wei - xlim[1])
      y <- x$pos[which(mins == min(mins))][1]
      return(y)
    }, xlim = xlim)
    heii <- lapply(my.inds2, function(x) {
      max(x[, cols])[1]
    })
    # saltiamo proposta picchi
    # if (suggested) {
    #   my.panel <- lapply(new.whole.data, function(popo) {
    #     pann <- big.peaks.col(popo$yy, tre = init.thresh)
    #     pann2 <- popo$xx[pann$pos]
    #     pann3 <- list(pos = pann$pos, hei = pann$hei,
    #                   wei = pann2)
    #     pkpn <- separate(pann3, type = "bp", shift = 1)
    #     return(list(wei = pkpn$wei, hei = pkpn$hei))
    #   })
    #   allpan <- unlist(lapply(my.panel, function(x) {
    #     x$wei
    #   }))
    #   allhei <- unlist(lapply(my.panel, function(x) {
    #     x$hei
    #   }))
    #   panel1.1 <- numeric()
    #   heis1.1 <- numeric()
    #   for (za in seq(1, 500, by = 1)) {
    #     step1 <- abs(za - allpan)
    #     good <- which(step1 < 0.48)
    #     if (length(good) > (length(n.inds) * 0.05)) {
    #       panel1.1[za] <- mean(allpan[good])
    #       heis1.1[za] <- mean(allhei[good])
    #     }
    #     else {
    #       panel1.1[za] <- NA
    #       heis1.1[za] <- NA
    #     }
    #   }
    #   if (is.null(xlim)) {
    #     panel.sugg <- panel1.1[-which(panel1.1 < min.panel |
    #                                     is.na(panel1.1))]
    #     heis.sugg <- heis1.1[-which(panel1.1 < min.panel |
    #                                   is.na(panel1.1))]
    #   }
    #   else {
    #     prov <- panel1.1[which(panel1.1 > xlim[1] & panel1.1 <
    #                              xlim[2])]
    #     bad <- which(is.na(prov))
    #     if (length(bad) > 0) {
    #       panel.sugg <- prov[-bad]
    #     }
    #     else {
    #       panel.sugg <- prov
    #     }
    #     prov2 <- heis1.1[which(panel1.1 > xlim[1] & panel1.1 <
    #                              xlim[2])]
    #     bad2 <- which(is.na(prov2))
    #     if (length(bad2) > 0) {
    #       heis.sugg <- prov[-bad2]
    #     }
    #     else {
    #       heis.sugg <- prov2
    #     }
    #   }
    # }
    tot.heii <- max(unlist(heii), na.rm = T)
    if (is.null(ylim)) {
      ylim <- c(0, tot.heii)
    }
    else {
      ylim <- ylim
    }
    layout(matrix(1, 1, 1))
    nn <- n.inds
    if (length(channel) > 1) {
      ylim[1] <- limosna[1]
      ylim[2] <- limosna[2]
    }
    plot(new.whole.data[[1]]$xx[-c(1:common[[1]])], y = new.whole.data[[1]]$yy[-c(1:common[[1]])], 
         type = "l", xlim = c(xlim[1], xlim[2]),  yaxt = "n", col = cols , xlab = "Size in base pairs", ylab = "DNA intensity in RFU", 
         xaxt = "n", lwd = lwd)
    #ylim = c(ylim[1], 
    #ylim[2])
    axis(1, at = seq(xlim[1], xlim[2], by = 2), labels = seq(xlim[1], 
                                                             xlim[2], by = 2), cex.axis = 0.7)
    #axis(2, at = seq(0, tot.heii, by = 500), labels = seq(0, 
    #                                                      tot.heii, by = 500), las = 1, cex.axis = 0.4)
    if (length(n.inds) == 1) {
      count <- count + 50
      setTxtProgressBar(pb, (count/tot))
    }
    else {
      count <- count + 1
    }
    if (length(n.inds) > 1) {
      for (i in 2:length(my.inds2)) {
        count <- count + 1
        a <- sum(unlist(heii)[1:(i - 1)])
        b <- sum(unlist(heii)[1:i])
        yy <- new.whole.data[[i]]$yy
        lines(new.whole.data[[i]]$xx[-c(1:common[[i]])], 
              y = yy[-c(1:common[[i]])], type = "l", col = transp(cfp[cols], 
                                                                  0.6), lwd = lwd)
        setTxtProgressBar(pb, (count/tot))
      }
    }
    # legend("topright", legend = "Peaks suggested", pch = 20, 
    #        col = "red", bty = "n", cex = 0.75)
    options(warn = oldw)
    close(pb)
  #   if (verbose) {
  #     cat("\n THE PEAKS RETURNED ARE SUGGESTIONS. \n   What you should do: \n a) Use the locator function, i.e. ''my.panel <- locator(type='p', pch=20, col='red')$x'' \n b) Click over the peaks you want to include in your panel \n c) Press the 'esc' key when done selecting peaks \n d) Make sure to provide the panel vector in the score.easy() function \n \n")
  #   }
  #   if (suggested) {
  #     points(x = panel.sugg, y = heis.sugg, pch = 20, cex = 0.7, 
  #            col = "red")
  #     points(x = panel.sugg, y = heis.sugg, cex = 0.9, 
  #            col = "black")
  #     suggested.list[[counter]] <- panel.sugg
  #   }
  #   par(new = TRUE)
  #   verbose = FALSE
  # }
  # par(new = FALSE)
  # if (suggested) {
  #   names(suggested.list) <- paste("channel_", channel, sep = "")
  #   return(suggested.list)
  # }
  layout(matrix(1, 1, 1))
  }
  return(new.whole.data) # right now 1 channel only! TODO implement list of all channels
}

overview2(my.inds=fd[1], channel = 1, ladder=ladder, lwd=2, channel.ladder=4, xlim=c(0, 1000), my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))

test <- getdata(my.inds=fd[1], channel = 1, ladder=ladder, lwd=2, channel.ladder=4, xlim=c(0, 1000), my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))
# sparizione segnale lì è per ylim, controllare quantili eblabla

pd2 <- data.frame(x=test[[1]]$xx, y=test[[1]]$yy)
ggplot(data=pd2, aes(y=y, x=x))+geom_line()+xlim(0,1000)

# e lo schifo prima del primo ladder?

write.table(pd2, '/tmp/test_channel1_fsa1.tsv', sep="\t", quote=F, row.names = F, col.names = T)

