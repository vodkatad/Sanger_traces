library(Fragman) 

fsa_file <- snakemake@input[['fsa']]
outplot_peak <- snakemake@output[['peaks']]
outplot_traces <- snakemake@output[['tracep']]
outf <- snakemake@output[['traces']]
ladder_chan <- as.numeric(snakemake@params[['ladder_chan']])
log_f <- snakemake@log[['log']]

wd <- getwd()
# From now on fragman code
# We create a temp dir with a single fsa since storing.inds loads all files from a given folder
tmp <- tempdir()
dir_name <- file.path(tmp, snakemake@wildcards[['sample']])
dir.create(dir_name)
file.copy(fsa_file, file.path(dir_name, basename(fsa_file)))
sink(log_f)
fd <- storing.inds(dir_name, lets.pullup=FALSE)
sink()
unlink(dir_name, recursive=TRUE)
setwd(wd)
# we change the default lets.pullup and set is as false
#lets.pullup
#A FALSE/TRUE value indicating if data should be treated for noise from channel to channel known as pull up or pull down peaks since wavelengths where the dyes are read usually overlap (blue->green->yellow->red->orange. The dafault is TRUE
# since we saw some of the previous signals to be 'strangely' modified 
#https://docs.google.com/presentation/d/1PMQWcPEmx1dVWgOFLy_crpC2k-jW654QFYbkU9bpXJ0/edit?usp=sharing
# slide 5-6-7
#in any case we shouldn't have fragments/peaks of the same lengths in the same mix

# Other steps: a smoothing Fourier transform and a correction for saturated peaks (code still not analyzed)

# These
# are the lengths of the ladder Fra used in the first experiment (I expect it to be the same but we need to check)
# right now hardcoded not an argument
ladder <-  unique(c(seq(60, 200, by=20), seq(200, 500, by=25), seq(500,600, by=50)))
# here we need to manually put channel.ladder=4 from the previous visual inspection
# trace has a trace:::find_ladders(d_fsas, config, show_progress_bar = TRUE) # questo aggiunge $trace_bp_df che è dove sono cercati i picchi.
# function but I'd start simple
## FIXME png and sink together?
png(outplot_peak, width =800, height = 400)
ladder.info.attach(stored=fd, ladder=ladder, channel.ladder=ladder_chan, prog=FALSE, draw=TRUE, attempt=15)
graphics.off()

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
#overview2(my.inds=fd[1], channel = 1:3, ladder=ladder, lwd=2, channel.ladder=4, xlim=c(80, 200), my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))
# in this case we have three automatically defined peaks only for channel_1
# parameters not all investigated
# launch this instruction and click on all the peaks you see, then press esc, this creates positions that will then be scored
#peaks <- locator(type='p', pch=20, col='red')$x
# locator is broken? it gets smaller x coordinates today, it was working when I worked on it previously
# resize of plot window in rstudio server? Yes, it gets easily broken, issues
# can be detected if red dots appear not where one clicked when using locator
#### mhhh.
#peaks <- c(99) # we could skip it and put the expected lengths of polyG on the genome/from Fra
# caution, it's really very sensitive (99 works, 100 doesn't)
#res <- score.markers(my.inds=fd[1], channel = 1, ladder=ladder, panel=peaks, electro=T, my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))

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
                     my.palette = NULL, verbose = TRUE, plotname= NULL) 
{
  png(plotname, width =800, height = 800)
  # controlliamo di avere un solo fsa, questo è forzato
  if (length(my.inds) != 1) {
    stop('We work with one fsa at a time')
  }
  ### indolore se non ci sono limiti li setta dai dati ####V
  if (length(channel) > 1 & is.null(ylim)) {
    limosna <- c(min(unlist(lapply(my.inds, min, na.rm = TRUE)), 
                     na.rm = TRUE), max(unlist(lapply(my.inds, max, na.rm = TRUE)), 
                                        na.rm = TRUE))
  }
  else {
    limosna <- ylim
  }
  res_df <- list()
  ############################################################V 
  suggested.list <- list()
  ############### itera su tutti i canali, ma in teoria noi gliene passiamo soltanto 1
  for (hhh in channel) {
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
    tot <- length(n.inds)
    ###############################V
    my.inds2 <- my.inds
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
    # list.models.inv <- lapply(list.data, function(da) {
    #   x <- da[[3]]
    #   y <- da[[1]]
    #   mod <- lm(y ~ x, data = da)
    #   return(mod)
    # })
    # i modelli girano sul canale del ladder ad ogni chiamata di overview per i diversi canali, vabbè.
    # xx per ogni oggetto fsa ha 1:sua lunghezza run, lavoriamo con 1 sola cols per volta dato il for su tutte a riga 118
    xx <- lapply(my.inds2, function(x, cols) {
      1:length(x[, cols])
    }, cols = cols)
    # xx è una lista di 1 elemento per noi che facciamo un fsa per volta
    # cols è il canale di questo giro
    newxx <- numeric()
    newyy <- numeric()
    new.whole.data <- list(NA)
    # Questo girerebbe per tutti i fsa, per noi no, fa un solo giro
    for (h in 1:length(xx)) {
      h1 <- n.inds[h]
      newxx <- as.vector(try(predict(list.models[[h1]], 
                                     newdata = data.frame(x = xx[[h]])), silent = TRUE))
      newyy <- my.inds2[[h]][, cols]
      ## questo è l'oggetto che vogliamo
      new.whole.data[[h]] <- list(xx = newxx, yy = newyy, origx=xx[[h]]) # aggiunta di origx
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
    tot.heii <- max(unlist(heii), na.rm = T)
    if (is.null(ylim)) {
      ylim <- c(0, tot.heii)
    }
    else {
      ylim <- ylim
    }
    nn <- n.inds
    if (length(channel) > 1) {
      ylim[1] <- limosna[1]
      ylim[2] <- limosna[2]
    }
    if (hhh==1) {
      print(plot(new.whole.data[[1]]$xx[-c(1:common[[1]])], y = new.whole.data[[1]]$yy[-c(1:common[[1]])], 
           type = "l", xlim = c(xlim[1], xlim[2]),  ylim = c(ylim[1], ylim[2]),
           yaxt = "n", col = transp(cfp[cols], 0.6) , xlab = "Size in base pairs", ylab = "DNA intensity in RFU", 
           xaxt = "n", lwd = lwd))
      print(axis(1, at = seq(xlim[1], xlim[2], by = 2), labels = seq(xlim[1], 
                                                                     xlim[2], by = 2), cex.axis = 0.7))
      print(axis(2, at = seq(0, tot.heii, by = 500), labels = seq(0, 
                                                                  tot.heii, by = 500), las = 1, cex.axis = 0.4))
      
    }  else {
      print(lines(new.whole.data[[1]]$xx[-c(1:common[[1]])], 
            y = new.whole.data[[1]]$yy[-c(1:common[[1]])],  type = "l", col = transp(cfp[cols], 
                                                                0.6), lwd = lwd))
    }
    options(warn = oldw)
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
    res_df[[hhh]] <- new.whole.data
  }
  graphics.off()
  return(res_df)
}

#overview2(my.inds=fd[1], channel = 1, ladder=ladder, lwd=2, channel.ladder=4, xlim=c(0, 1000), my.palette=c('#FFFF00', '#00F5FF','#7FFF00'))

# right only if ladder is the last channel FIXME
lastch <- ladder_chan
list_of_list <- getdata(my.inds=fd, channel = 1:lastch, ladder=ladder, lwd=2, channel.ladder=ladder_chan,
                        my.palette=c('#FFFF00', '#00F5FF','#7FFF00','red'), plotname=outplot_traces, verbose=FALSE)


# unchecked assumptions FIXME
res <- list()
for (ch in 1:lastch) {
  d <- list_of_list[[ch]]
  if (ch == 1) {
    res[['mw']] <- d[[1]]$xx
    res[['orig_x']] <- d[[1]]$origx
  }
  res[[paste0('channel_', ch)]] <- d[[1]]$yy
  #pd2 <- data.frame(mw=d[[1]]$xx, intensity=d[[1]]$yy, orig_x=d[[1]]$origx)
}
write.table(as.data.frame(res), outf, sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

