
#' @title Setup the watchdog system
#' 
#' @description This function generates a number of bach scripts and a template config file for the watchdog
#' 
#' @param bindir The bash scripts are written in this directory
#' @param hdf5dir The directory on the server for writing the HDF5 files containing the full images
#' @param hdf5validationdir The directory on the server for writing the status files indicating if the TIF images and the HDF5 files are identical
#' @param hdf5projection The directory on the server for writing the HDF5 files containing the z-projection
#' @param htmldir The directory on the server for writing thumbnail images
#' @param layoutdir The directory on the server for the plate layouts
#' @param featuresdir The directory on the server for the feature files
#' @param configdir The configuration file, which defines all the directories, the number of z-stacks, channels, and fields
#' @param segmentationdir The directory in which the DNN segmentation is stored
#' @param segmentationhtmldir The directory in which the preview images for the segmentation are stored
#' @param pythondir The directory that contains the python script(s) for the segmentation as well as the weights and layout information for dnnSwift
#' @param tifdiagnosticsdir The directory in which TIF diagnostic files are stored
#' @param showWarnings A boolean value indicating if a warning should be shown if the config file already exists when this function is executed
#' 
#' @return Nothing is returned
#' 
#' @author Bernd Fischer
#' 
#' @examples print(setupWatchdog)
#' @export
setupWatchdog <- function(bindir = "~/bin",
                          pythondir = file.path(getwd(), "config/segmentation_config"),
                          hdf5dir = file.path(getwd(),"hdf5dir"), 
                          hdf5validationdir = file.path(getwd(),"hdf5validation"),
                          hdf5projection = file.path(getwd(),"hdf5projection"), 
                          htmldir = file.path(getwd(),"htmldir"),
                          layoutdir = file.path(getwd(),"layouts"),
                          featuresdir = file.path(getwd(), "features"),
                          segmentationdir = file.path(getwd(), "segmentation"),
                          segmentationhtmldir = file.path(getwd(), "segmentationhtmldir"),
                          tifdiagnosticsdir = file.path(getwd(), "tifdiagnostics"),
                          configdir = file.path(getwd(),"configdir"),
                          showWarnings = TRUE) {
    dir.create(configdir, showWarnings = FALSE)
    if (!file.exists(file.path(configdir, "watchdogConfig.R"))) {
        script = c("sleepSeconds=300",
                   "nrWells=384",
                   sprintf("indir=<__SET MANUALLY__>"),
                   sprintf("pythondir=\"%s\"", pythondir),
                   sprintf("hdf5dir=\"%s\"",hdf5dir),
                   sprintf("hdf5validationdir=\"%s\"",hdf5validationdir),
                   sprintf("hdf5projection=\"%s\"",hdf5projection),
                   sprintf("htmldir=\"%s\"",htmldir),
                   sprintf("layoutdir=\"%s\"",layoutdir),
                   sprintf("featuresdir=\"%s\"",featuresdir),
                   sprintf("segmentationdir=\"%s\"",segmentationdir),
                   sprintf("tifdiagnosticsdir=\"%s\"",tifdiagnosticsdir),
                   sprintf("segmentationhtmldir=\"%s\"",segmentationhtmldir),
                   "pythonexec=\"PROMISE_DNNSegmentation.py\"",
                   "qsub=\"qsub -q ag_fischer -l walltime=18:00:00,select=1:ncpus=1:mem=10GB -W umask=0011\"",
                   "startScript=\"library(PROMISE);startWorkflow();warnings()\"",
                   "fctcall=\"PROMISEworkflow\"",
                   "channels=c('Cy3','FITC','DAPI')",
                   "exciteWavelengths=c('Green', 'Blue', 'UV')",
                   "fields=c('1','2','3','4','5','6','7','8','9')",
                   "fields_layout=c('1','2','3','6','5','4','7','8','9')",
                   "stacks = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11')",
                   "hdf5_compression=3",
                   "hdf5_chunksize=256")
        writeLines(script, file.path(configdir, "watchdogConfig.R"))
        dir.create(hdf5dir, showWarnings = FALSE, recursive = TRUE)
        dir.create(hdf5projection, showWarnings = FALSE, recursive = TRUE)
        dir.create(htmldir, showWarnings = FALSE, recursive = TRUE)
        dir.create(featuresdir, showWarnings = FALSE, recursive = TRUE)
        dir.create(hdf5validationdir, showWarnings = FALSE, recursive = TRUE)
        dir.create(segmentationdir, showWarnings = FALSE, recursive = TRUE)
        dir.create(tifdiagnosticsdir, showWarnings = FALSE, recursive = TRUE)
        dir.create(segmentationhtmldir, showWarnings = FALSE, recursive = TRUE)
        message("Write config file ",file.path(configdir, "watchdogConfig.R"))
    } else {
        if (showWarnings) {
            warning("config file already exists.")
        }
    }

    dir.create(bindir, showWarnings = FALSE)
    writeLines(c("#!/bin/bash",
                 "",
                 sprintf('R -e \'library(PROMISE); watchdogStart(configdir=\"%s\")\' &',
                 configdir),
                 "echo \"\"",
                 ""),file.path(bindir,"watchdogStart.sh"))
    system(sprintf("chmod a+x %s",file.path(bindir,"watchdogStart.sh")))
    message("Write start script ",file.path(bindir,"watchdogStart.sh"))

    writeLines(c("#!/bin/bash",
                 "",
                 sprintf('R -e \'library(PROMISE); watchdogStop(configdir=\"%s\")\'',
                         configdir),
                 "echo \"\"",
                 ""),file.path(bindir,"watchdogStop.sh"))
    system(sprintf("chmod a+x %s",file.path(bindir,"watchdogStop.sh")))
    message("Write stop script ",file.path(bindir,"watchdogStop.sh"))

    writeLines(c("#!/bin/bash",
                 "",
                 "echo ''",
                 "echo 'last loop of watchdog started at'",
                 sprintf('tail -n 20 %s',file.path(configdir,"status","lastloopstart.txt")),
                 "echo ''",
                 sprintf('echo "last 20 lines of log file %s"',file.path(configdir,"status","log.txt")),
                 "echo ''",
                 sprintf('tail -n 20 %s',file.path(configdir,"status","log.txt"))),
               file.path(bindir,"watchdogLog.sh"))
    system(sprintf("chmod a+x %s",file.path(bindir,"watchdogLog.sh")))
    message("Write log script ",file.path(bindir,"watchdogLog.sh"))
    
    dir.create(bindir, showWarnings = FALSE)
    writeLines(c("#!/bin/bash",
                 "",
                 sprintf('R -e \'library(PROMISE); processWell()\' --args $1 $2 $3 %s &',configdir),
                 "echo \"\"",
                 ""),file.path(bindir,"processWell.sh"))
    system(sprintf("chmod a+x %s",file.path(bindir,"processWell.sh")))
    message("Write processWell script ",file.path(bindir,"processWell.sh"))

    invisible(NULL)
}


#' @title Stop the watchdog
#' 
#' @description This function kills the process of the R session running the watchdog. You have to call watchdogStop before restarting it after a crash.
#' 
#' @param configdir The configuration file, which defines all the directories, the number of z-stacks, channels, and fields
#' 
#' @return Nothing is returned
#' 
#' @author Bernd Fischer
#' 
#' @examples print(watchdogStop)
#' @export
watchdogStop <- function(configdir) {
    if (file.exists(file.path(configdir, "status","pid.txt"))) {
        logfile = file(file.path(configdir, "status","log.txt"),open = "a")
        writeLines(paste(Sys.time(),"Stop watching directory"),con=logfile)
        close(logfile)
        PID = readLines(file.path(configdir, "status","pid.txt"))
        system(paste0("kill -9 ",PID))
        message("process with PID ",PID, " killed.\n")
        file.remove(file.path(configdir, "status","pid.txt"))
    } else {
        message("No running watchdog detected.")
    }
    invisible(NULL)
}

#' @title Start the watchdog
#' 
#' @description Starts the watchdog. Before starting the watchdog the first time, you have to set it up with setupWatchdog. Note that you have to call watchdogStop before restarting it after a crash.
#' 
#' @param configdir The configuration file, which defines all the directories, the number of z-stacks, channels, and fields
#' 
#' @return Nothing is returned
#' 
#' @author Bernd Fischer
#' 
#' @examples print(watchdogStart)
#' @export
watchdogStart <- function(configdir) {
    if (!file.exists(file.path(configdir, "watchdogConfig.R"))) {
        stop("file watchdogConfig.R not found in directory ",configdir,". Run setupWatchdog() to generate one")
    } else {
        source(file.path(configdir, "watchdogConfig.R"))
    }

    dir.create(file.path(configdir, "status"), recursive = TRUE, showWarnings = TRUE)
    dir.create(file.path(configdir, "status", "processedWells"), recursive = TRUE, showWarnings = TRUE)

    # check if another watchdog is running
    if (file.exists(file.path(configdir, "status","pid.txt"))) {
        stop("It is likely that another watchdog is still running. Call watchdogStop() to stop it.")
    }
    writeLines(as.character(Sys.getpid()),
               file.path(configdir, "status","pid.txt"))

    # if (!file.exists(file.path(configdir, "status","processedPlates.txt"))) {
    file.create(file.path(configdir, "status","processedPlates.txt"), showWarnings = FALSE) # Requeue all plates in indir
    # }
    if (!file.exists(file.path(configdir, "status","stopProcessingPlates.txt"))) {
        file.create(file.path(configdir, "status","stopProcessingPlates.txt"), showWarnings = FALSE)
    }
    if (!file.exists(file.path(configdir, "status","log.txt"))) {
        file.create(file.path(configdir, "status","log.txt"), showWarnings = FALSE)
    }
    if (!file.exists(file.path(configdir, "status","lastloopstart.txt"))) {
        file.create(file.path(configdir, "status","lastloopstart.txt"), showWarnings = FALSE)
    }
    if (!file.exists(file.path(configdir, "status","cmd.R"))) {
        file.create(file.path(configdir, "status","cmd.R"), showWarnings = FALSE)
    }
    AllNoBarcode = c()
    logfile = file(file.path(configdir, "status","log.txt"),open = "a")
    writeLines(paste(Sys.time(),"Start watching directory ",indir),con=logfile)
    flush(logfile)
    W = wells(nrWells)
    queue = data.frame(platename = c(), plateIndir = c(), well=c(),stringsAsFactors = FALSE)
    while (TRUE) {
        writeLines(paste(Sys.time()," (start of loop)"),con=file.path(configdir, "status","lastloopstart.txt"))
        if (file.exists(file.path(configdir, "cmd.R"))) {
            source(file.path(configdir, "cmd.R"))
        }
        if (file.exists(file.path(configdir, "cmdOnce.R"))) {
            e = try({ source(file.path(configdir, "cmdOnce.R")) })
            file.remove(file.path(configdir, "cmdOnce.R"))
        }

        NPBS = length(grep("ag_fischer",system("qstat", intern = TRUE)))
        if (NPBS <= 200) {
            # check for new plates
        processedPlates = readLines(file.path(configdir, "status","processedPlates.txt"))
        platesIndir = dir(path = indir)
        hasBarcode = !grepl("Plate",platesIndir)
        # Extract plate barcode from directory name
        gr = gregexpr(pattern = "_",text = platesIndir)
        hasBarcode = hasBarcode & (sapply(gr, length) == 2)
        NoBarcode = platesIndir[!hasBarcode]
        platesIndir = platesIndir[hasBarcode]
        gr = gregexpr(pattern = "_",text = platesIndir)
        platenames = platesIndir
        for (i in seq_along(platesIndir)) {
            platenames[i] = substr(platesIndir[i],gr[[i]][1]+1,gr[[i]][2]-1)
        }
        I = which(!(platenames %in% processedPlates))
        platenames = platenames[I]
        platesIndir = platesIndir[I]

        # check if barcode is read properly, Do not process plates with platenames containing "Plate"
        NoBarcode = NoBarcode[!(NoBarcode %in% AllNoBarcode)]
        if (length(NoBarcode)) {
            AllNoBarcode = c(AllNoBarcode,NoBarcode)
            writeLines(paste(Sys.time(),"Plate ",NoBarcode," without barcode detected. Plate will not be processed."),con=logfile)
            flush(logfile)
        }

        n = length(platesIndir)
        if (n > 0) {
            for (i in seq_along(platesIndir)) {
                p = platenames[i]
                dfnew = data.frame(platename = platenames[i], plateIndir = platesIndir[i], 
                                   row = W$rows, col = W$cols, stringsAsFactors = FALSE)
                # Do not add wells to queue that are already send to qsub
                if (file.exists(file.path(configdir, "status","processedWells",p))) {
                    processedWells = readLines(file.path(configdir, "status","processedWells",p))
                    dfnew = dfnew[!(paste(dfnew$platename, dfnew$row, dfnew$col,sep="_") %in% 
                                        processedWells),,drop=FALSE]
                }
                if (!file.exists(file.path(configdir, "status","processedWells",p))) {
                    file.create(file.path(configdir, "status","processedWells",p), showWarnings = FALSE)
                }

                # Remove all entries with plate p from watching queue
                queue = queue[which(queue$platename != p),,drop=FALSE]

                queue = rbind(queue,dfnew)
                f = file(file.path(configdir, "status","processedPlates.txt"),open="a")
                writeLines(p, f)
                close(f)
                writeLines(paste(Sys.time(),"Start processing plate",p),con=logfile)
                
                # Run the rename. This is here instead of in a job script as it needs
                # to read many files from the tiff-directory. This could cause problems
                # if it happens in parallel. Plus the function only generates symlinks
                # so it's low-latency to begin with.
                renameTIFF(plateIndir = platesIndir[i], configdir = configdir)
                
                flush(logfile)
            }
        }

        # Stop processing of plates 
        stopProcessingPlates = readLines(file.path(configdir, "status","stopProcessingPlates.txt"))
        if(any(queue$platename %in% stopProcessingPlates)) {
            writeLines(paste(Sys.time(),"Stop processing plate ",queue$platename[queue$platename %in% stopProcessingPlates]),con=logfile)
            queue = queue[!(queue$platename %in% stopProcessingPlates),,drop=FALSE]
            flush(logfile)
        }

        # Do not process wells that already have been computed
        # possibleHdf5Files = apply(queue, 1, function(x) hdf5Filename(hdf5dir, x["platename"], x["row"], x["col"], configdir))
        # I = which(file.exists(possibleHdf5Files))
        # if (length(I) > 0) {
        #     writeLines(paste(Sys.time(),"Do not process plate ",queue$platename[I],
        #                      "well",queue$row[I],queue$col[I],"file already exists"),
        #                con=logfile)
        #     flush(logfile)
        #     queue = queue[-I,,drop=FALSE]
        # }

#        save(queue, file=file.path(configdir,"status","queue.rda"))

        if (nrow(queue) > 0) {
            P = unique(queue$platename)
            for (p in P) {
                dir.create(file.path(configdir,"status","log",p), showWarnings = FALSE, recursive = TRUE)
            }
        }

        # Start processing wells
        removeFromQueue = c()
        for (i in seq_len(nrow(queue))) {
            if (NPBS <= 250) {
            if (filesComplete(plateIndir=queue$plateIndir[i], 
                              platename=queue$platename[i], 
                              row=queue$row[i], col=queue$col[i],configdir = configdir)) {
                submitClusterQueue(queue$plateIndir[i], queue$platename[i], queue$row[i], queue$col[i],
                                   configdir=configdir, logfile=logfile)
                removeFromQueue = c(removeFromQueue,i)
            }
            NPBS = NPBS + 1
            }
        }
        if (length(removeFromQueue) > 0) {
            queue = queue[-removeFromQueue,,drop=FALSE]
        }
        } # check for NPBS
        Sys.sleep(sleepSeconds)
    }
}

#' @title Submits a job to the cluster queue to process one well
#' 
#' @description This function submits a job to the cluster queue to process one well. 
#' 
#' @param plateIndir The base directory on the server, which contains a folder for each plate
#' @param platename The ID (and folder name) of the plate
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The configuration file, which defines the number of z-stacks, channels, and fields
#' @param logfile An open connection to the logfile. If missing the logfile is opened and closed again after adding a log entry
#' @param custom_function An optional custom workflow command. Must be one of the PROMISE* functions.
#' 
#' @return A boolean value indicating if the projection was successful
#' 
#' @author Bernd Fischer
#' 
#' @examples print(submitClusterQueue)
#' @export
submitClusterQueue <- function(plateIndir, platename, row, col, configdir, logfile, custom_function = NULL) {
    if (!file.exists(file.path(configdir, "watchdogConfig.R"))) {
        stop("file watchdogConfig.R not found in directory ",configdir,". Run setupWatchdog() to generate one")
    } else {
        source(file.path(configdir, "watchdogConfig.R"))
    }
    if (missing(plateIndir)) {
        D = dir(indir)
        I = grep(platename,D)
        if (length(I) == 0) {
            warning("incoming directory not found. If incoming directory is expected for the function call, the process will fail.")
            plateIndir = "unknownDirectory"
        } else {
            if (length(I) > 1) {
                stop("multiple incoming directories found for this platename")
            } else {
                plateIndir = D[I]
            }
        }
    }
    # Inject custom function call
    if(!is.null(custom_function)) {
      workflows = grep("^PROMISE", as.vector(lsf.str("package:PROMISE")), value = TRUE)
      if(custom_function %in% workflows) {
        fctcall = custom_function
      }
    }
    file.remove(file.path(configdir, "status", "log", platename,
                          sprintf("%s_%s_%s_%s_err.txt",platename,row, col,fctcall)))
    file.remove(file.path(configdir, "status", "log", platename,
                          sprintf("%s_%s_%s_%s_out.txt",platename,row, col,fctcall)))
    # shellscript = sprintf("echo \'export R_LIBS_USER=\"/collab-ag-fischer/Rlibs\"; R --vanilla -e \"%s\" --args %s %s %s %s %s %s\' | %s -e %s -o %s",
    shellscript = sprintf("echo \'R --vanilla -e \"%s\" --args %s %s %s %s %s %s\' | %s -e %s -o %s",
                          startScript,
                          fctcall, plateIndir,
                          platename,
                          row, col,
                          configdir,
                          qsub,
                          file.path(configdir, "status", "log", platename,
                                    sprintf("%s_%s_%s_%s_err.txt",platename,row, col,fctcall)),
                          file.path(configdir, "status", "log", platename,
                                    sprintf("%s_%s_%s_%s_out.txt",platename,row, col,fctcall)))
    res = system(shellscript, intern = TRUE)
    f = file(file.path(configdir, "status","processedWells",platename),open="a")
    writeLines(paste(platename,row,col,collapse="_",sep="_"),f)
    close(f)
    if (missing(logfile)) {
        logfile = file(file.path(configdir, "status","log.txt"),open = "a")
        writeLines(paste(Sys.time(),"Restart processing plate",platename,"well",
                         row, col),con=logfile)
        writeLines(shellscript,con=logfile)
        writeLines(res,con=logfile)
        close(logfile)
    } else {
        writeLines(paste(Sys.time(),"Start processing plate",platename,"well",
                         row, col),con=logfile)
        writeLines(shellscript,con=logfile)
        writeLines(res,con=logfile)
        flush(logfile)
    }
    return(TRUE)
}


#' @title Read command line and start R function
#' 
#' @description This function is called by the watchdog to start the processing
#'  of a single well. It reads the arguments from command line and passes them
#'  to the function specified as its first argument. This function is only 
#'  usefull, if you want to process a single well by a shell command. You can
#'  do so by calling
#'  >R -e "library(PROMISE);startWorkflow --args tif2hdf5 plateIndir platename row col configdir".
#' 
#' @return Nothing is returned
#' 
#' @author Bernd Fischer
#' 
#' @examples print(startWorkflow)
#' @export
startWorkflow <- function() {
    cmdArgs = commandArgs(trailingOnly = TRUE)
    do.call(cmdArgs[1],list(plateIndir=cmdArgs[2], 
                            platename=cmdArgs[3],
                            row=cmdArgs[4],
                            col=cmdArgs[5],
                            configdir=cmdArgs[6]))
    invisible(NULL)
}



#' @title Submits a job to the cluster queue to process a single well
#' 
#' @description This function is called by the batch script processWell.sh.
#'  It submits a single well process to the cluster queue.
#' 
#' @return Nothing is returned
#' 
#' @author Bernd Fischer
#' 
#' @examples print(processWell)
#' @export
processWell <- function() {
    cmdArgs = commandArgs(trailingOnly = TRUE)
    do.call("submitClusterQueue",list(platename=cmdArgs[1],
                            row=cmdArgs[2],
                            col=cmdArgs[3],
                            configdir=cmdArgs[4]))
    invisible(NULL)
}

#' @title Test if all files of one well are complete
#' 
#' @description Tests if all required files for processing one well are 
#' existing and not open any more.
#' 
#' @param plateIndir The base directory on the server, which contains a folder for each plate
#' @param platename The ID (and folder name) of the plate
#' @param row The row of the well
#' @param col The column of the well
#' @param configdir The configuration file, which defines all the directories, the number of z-stacks, channels, and fields
#'
#' @return Logical indicating if all files are complete
#' 
#' @author Bernd Fischer
#' 
#' @examples print(filesComplete)
#' @export
filesComplete <- function(plateIndir, platename, row, col, configdir) {
    fn <- imageFilenames(plateIndir, platename, row, col, configdir)
    res = FALSE
    if (all(file.exists(fn))) {
        if (all(!file.isopen(fn))) {
            res = TRUE
        }
    }
    res
}

#' @title Test if files are open
#' 
#' @description Tests for each file, if it is open
#' 
#' @param fn A vector of filenames
#'
#' @return Logical indicating for each file, if it is open
#' 
#' @author Bernd Fischer
#' 
#' @examples print(file.isopen)
#' @export
file.isopen <- function(fn) {
    res = rep(FALSE, length(fn))
    for (i in seq_along(fn)) {
        if (length(suppressWarnings(system(sprintf("lsof -- '%s'",fn[i]), intern = TRUE)) > 0)) {
            res[i] = TRUE
        }
    }
    res
}

# indir = "/Users/fischebe/Documents/SVN/fischerlab/users/sauerja/Projects/PROMISE/PROMISE/vignette/indir"
# hdf5dir = "/Users/fischebe/Documents/SVN/fischerlab/users/sauerja/Projects/PROMISE/PROMISE/vignette/outdir"
# configdir = "/Users/fischebe/Documents/SVN/fischerlab/users/sauerja/Projects/PROMISE/PROMISE/vignette/config"
