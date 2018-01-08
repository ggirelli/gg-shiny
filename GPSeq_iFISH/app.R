
# --------------------------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: generate plots for GPSeq and iFISH comparison.
# 
# --------------------------------------------------------------------------------------------------



# DEPENDENCIES =====================================================================================
library(shiny)
library(plotly)
library(gridExtra)
library(ggplot2)

# FUNCTIONS ========================================================================================

order_allele_labels = function(data) {
    data = factor(data, levels = c("Extra", "Single", "Central", "Peripheral"))
    return(data)
}

order_chromosomes = function(data) {
    chrid = unlist(lapply(levels(data), FUN = function(x) substr(x, 4, nchar(x)) ))
    chrid[chrid == "X"] = 23
    chrid[chrid == "X"] = 24
    chrid = as.numeric(chrid)
    return(factor(data, levels = levels(data)[order(chrid)]))
}

order_data = function(data, flag) {
    if ( flag == "Allele" ) { return(order_allele_labels(data)) }
    if ( flag == "chr" ) { return(order_chromosomes(data)) }
    return(data)
}

filter_cell_type = function(data, input) {
    # FIlters data based on cell_type column and cell_type selected in input$cell_type
    # If input$all_cell_types checkbox is selected, this is skipped.
    if( input$all_cell_types ) return(data)
    return(data[data$cell_type == input$cell_type,])
}

filter_dataset = function(data, input) {
    # FIlters data based on dataset column and dataset selected in input$dataset_label.
    # If input$all_datasets checkbox is selected, this is skipped.
    dataset = unlist(strsplit(input$dataset_label, '/', fixed = T))[1]
    if( input$all_datasets ) return(data)
    return(data[data$dataset == dataset,])
}

filter_cell_phase = function(data, input, flag) {
    # Filters data based on cell cycle phase and input$cell_phase selection.
    # For G1: based on G1 column.
    # For 1 allele:
    # For 2 alleles:

    data$uniID = paste0(
        data$dataset, "~",
        data$File, "~",
        data$cell_ID
    )

    # Select 1-allele cells
    if( "Only cells with 1 allele" == input$cell_phase ) {
        if( "alleles" == flag ) return()
        if( "nuclei" != flag ) {
            fpath = paste0("data/dots.only1allele.tsv")
            if ( !file.exists(fpath) ) {
                data = data[0 != data$cell_ID,]
                data = do.call(rbind, by(data, data$uniID,
                    FUN = function(subt) {
                        if( all(1 == table(subt$Channel)) ) return(subt)
                    }
                ))

                write.table(data, fpath, quote = F, row.names = F, col.names = T, sep = "\t")
            } else {
                data = read.delim(fpath, as.is = T, header = T)
            }
        }
    }

    # Select 2-alleles cells
    fpath = paste0("data/dots.only2allele.tsv")
        if( ! flag %in% c("alleles", "nuclei") ) {
            if( "Only cells with 2 alleles" == input$cell_phase ) {
            if ( !file.exists(fpath) ) {
                data = data[0 != data$cell_ID,]
                data = do.call(rbind, by(data, data$uniID,
                    FUN = function(subt) {
                        if( all(2 == table(subt$Channel)) ) return(subt)
                    }
                ))

                write.table(data, fpath, quote = F, row.names = F, col.names = T, sep = "\t")
            } else {
                data = read.delim(fpath, as.is = T, header = T)
            }
        }
    }
    
    # Select only G1
    if( T == input$g1only ) return(data[data$G1 == 1,])

    # No selection
    return(data)
}

select_dot_data = function(input) {
    # Applies G1 and dataset selection to single-dot data.
    data = filter_cell_phase(md, input,  "dots")
    data = filter_cell_type(data, input)
    data = filter_dataset(data, input)

    return(data)
}

select_allele_data = function(input) {
    # Applies G1 and dataset selection to allele data.
    data = filter_cell_phase(mda, input, "alleles")
    data = filter_cell_type(data, input)
    data = filter_dataset(data, input)
    return(data)
}

select_nuclei_data = function(input) {
    # Applies dataset selection to allele data.
    data = filter_cell_phase(nd, input, "nuclei")
    data = filter_cell_type(data, input)
    data = filter_dataset(data, input)
    return(data)
}

label2col = function(label, set) {
    return(names(which(set == label))[1])
}

# PREPARE DATA =====================================================================================

# Read main data
dataset_date = "2018-01-03"
md = read.delim("data/2018-01-03_dots.merged.tsv", as.is = T, header = T)
mda = read.delim("data/2018-01-03_alleles.merged.tsv", as.is = T, header = T)
nd = read.delim("data/2018-01-03_nuclei.merged.tsv", as.is = T, header = T)
cs = read.delim("data/chr_size.txt", as.is = T, header = F)
colnames(cs) = c("chr", "size")

# Convert allele labels to readable ones
md$Allele[md$Allele == -1] = "Extra"
md$Allele[md$Allele == 0] = "Single"
md$Allele[md$Allele == 1] = "Central"
md$Allele[md$Allele == 2] = "Peripheral"

# Add chromosome, size and set columns
md$chr = unlist(lapply(strsplit(md$label, '.', fixed = T), FUN = function(x) x[[1]]))
md$chr_size = cs$size[match(tolower(md$chr), cs$chr)]
md$set = unlist(lapply(strsplit(md$label, '.', fixed = T), FUN = function(x) x[[2]]))
mda$chr = unlist(lapply(strsplit(mda$label, '.', fixed = T), FUN = function(x) x[[1]]))
mda$chr_size = cs$size[match(tolower(mda$chr), cs$chr)]
mda$set = unlist(lapply(strsplit(mda$label, '.', fixed = T), FUN = function(x) x[[2]]))

# Make channel lower-case
md$Channel = tolower(md$Channel)
mda$Channel = tolower(mda$Channel)

# Add probe name columns
md$probe_name = paste0(md$label, '.', md$Channel)

# Column type
md_col_type = list(
    "File" = "factor",
    "Channel" = "factor",
    "Nuclei" = "factor",
    "x" = "real",
    "y" = "real",
    "z" = "real",
    "Value" = "skip",
    "FWHM" = "skip",
    "Label" = "factor",
    "cell_ID" = "factor",
    "lamin_dist" = "real",
    "lamin_dist_norm" = "real",
    "centr_dist" = "real",
    "centr_dist_norm" = "real",
    "dilation" = "skip",
    "angle" = "skip",
    "com" = "skip",
    "version" = "skip",
    "G1" = "factor",
    "Allele" = "factor",
    "dataset" = "factor",
    "cell_type" = "factor",
    "label" = "factor",
    "chr" = "factor",
    "chr_size" = "real",
    "set" = "factor",
    "probe_name" = "factor"
)
mda_col_type = list(
    "File" = "factor",
    "Channel" = "factor",
    "cell_ID" = "factor",
    "G1" = "factor",
    "d_3d" = "real",
    "d_lamin" = "real",
    "d_lamin_norm" = "real",
    "d_centr" = "real",
    "d_centr_norm" = "real",
    "angle" = "real",
    "dataset" = "factor",
    "label" = "factor",
    "cell_type" = "factor",
    "chr" = "factor",
    "chr_size" = "real",
    "set" = "factor"
)
nd_col_type = list(
    "s" = "factor",
    "n" = "factor",
    "flat_size" = "real",
    "size" = "real",
    "surf" = "real",
    "sumI" = "real",
    "meanI" = "real",
    "shape" = "real",
    "G1" = "factor",
    "dataset" = "factor",
    "cell_type" = "factor"
)

# Column label
md_col_label = list(
    "File" = "Field",
    "Channel" = "Channel",
    "Nuclei" = "Cell ID (by DOTTER)",
    "x" = "X-coord",
    "y" = "Y-coord",
    "z" = "Z-coord",
    "Value" = NA,
    "FWHM" = NA,
    "Label" = "Allele label (by DOTTER)",
    "cell_ID" = "Cell ID (by GPSeq)",
    "lamin_dist" = "Absolute distance from lamina [nm]",
    "lamin_dist_norm" = "Relative distance from lamina [a.u.]",
    "centr_dist" = "Absolute distance from central region [nm]",
    "centr_dist_norm" = "Relative distance from central region [a.u.]",
    "dilation" = NA,
    "angle" = NA,
    "com" = NA,
    "version" = NA,
    "G1" = "G1-status",
    "Allele" = "Allele label (by GPSeq)",
    "dataset" = "Dataset ID",
    "cell_type" = "Cell type",
    "label" = "Probe set name",
    "chr" = "Chromosome",
    "chr_size" = "Chromosome length [bp]",
    "set" = "Probe set type",
    "probe_name" = "Probe name"
)
mda_col_label = list(
    "File" = "Field",
    "Channel" = "Channel",
    "cell_ID" = "Cell ID (by GPSeq)",
    "G1" = "G1-status",
    "d_3d" = "Absolute 3D distance [nm]",
    "d_lamin" = "Distance between absolute lamina layers [nm]",
    "d_lamin_norm" = "Distance between normalized lamina layers [a.u.]",
    "d_centr" = "Distance between absolute center layers [nm]",
    "d_centr_norm" = "Distance between normalized center layers [a.u.]",
    "angle" = "Angle with nucleus center of mass [deg]",
    "dataset" = "Dataset ID",
    "label" = "Probe set name",
    "cell_type" = "Cell type",
    "chr" = "Chromosome",
    "chr_size" = "Chromosome length [bp]",
    "set" = "Probe set type"
)
nd_col_label = list(
    "s" = "Field",
    "n" = "Cell ID (by GPSeq)",
    "flat_size" = "Nuclear flattened area [px]",
    "size" = "Nuclear volume [vx]",
    "surf" = "Nuclear surface [px]",
    "sumI" = "Intensity integral over nuclear volume [a.u.]",
    "meanI" = "Average intensity over nuclear volume [a.u.]",
    "shape" = "Sphericity",
    "G1" = "G1-status",
    "dataset" = "Dataset ID",
    "cell_type" = "Cell type"
)

# Column description
md_col_descr = list(
    "File" = "ID of the field of view.",
    "Channel" = "Label of a channel.",
    "Nuclei" = "Nucleus ID assigned by DOTTER.",
    "x" = "X-coordinate of a single FISH signal, by DOTTER.",
    "y" = "Y-coordinate of a single FISH signal, by DOTTER.",
    "z" = "Z-coordinate of a single FISH signal, by DOTTER.",
    "Label" = "Allele label assigned by DOTTER.",
    "cell_ID" = "Nucleus ID assigned by gpseq-img-py.",
    "lamin_dist" = "Absolute distance from nuclear lamina [nm].",
    "lamin_dist_norm" = "Relative distance from nuclear lamina [a.u.].",
    "centr_dist" = "Absolute distance from nucleus central region [nm].",
    "centr_dist_norm" = "Relative distance from nucleus central region [a.u.].",
    "G1" = "1 for G1 cells, 0 for non-G1 cells.",
    "Allele" = "Allele labeling by gpseq-img-py.
    Extra: more than 2 alleles in that cell.
    Single: only one allele in that cell.
    Central: central allele, from a cell with 2 alleles.
    Peripheral: peripheral allele, from a cell with 2 alleles.",
    "dataset" = "Dataset ID.",
    "cell_type" = "Cell type.",
    "label" = "Probe set name.",
    "chr" = "Chromosome.",
    "chr_size" = "Chromosome length in bp.",
    "set" = "Probe set type (c: central; s: start; e: end).",
    "probe_name" = "Probe name (chr.type.channel)."
)
mda_col_descr = list(
    "File" = "ID of the field of view.",
    "Channel" = "Label of a channel.",
    "cell_ID" = "Nucleus ID assigned by gpseq-img-py.",
    "G1" = "1 for G1 cells, 0 for non-G1 cells.",
    "d_3d" = "3D distance between two alleles in the same channel.",
    "d_lamin" = "Distance between absolute lamin layers of two alleles in the same channel.",
    "d_lamin_norm" = "Distance between relative lamin layers of two alleles in the same channel.",
    "d_centr" = "Distance between absolute center layers of two alleles in the same channel.",
    "d_centr_norm" = "Distance between relative center layers of two alleles in the same channel.",
    "angle" = "Angle between two alleles in the same channel and the nucleus center of mass",
    "dataset" = "Dataset ID.",
    "label" = "Probe set name.",
    "cell_type" = "Cell type.",
    "chr" = "Chromosome.",
    "chr_size" = "Chromosome length in bp.",
    "set" = "Probe set type (c: central; s: start; e: end)."
)
nd_col_descr = list(
    "s" = "ID of the field of view.",
    "n" = "Nucleus ID assigned by gpseq-img-py.",
    "flat_size" = "Nuclear flattened area [px]",
    "size" = "Nuclear volume [vx]",
    "surf" = "Nuclear surface [px]",
    "sumI" = "Intensity integral over nuclear volume [a.u.]",
    "meanI" = "Average intensity over nuclear volume [a.u.]",
    "shape" = "Sphericity as ratio between the surface of the nucleus and of a sphere with the same
        volume.",
    "G1" = "1 for G1 cells, 0 for non-G1 cells.",
    "dataset" = "Dataset ID.",
    "cell_type" = "Cell type."
)

# Summary type to function
sum_types = list(
    "Variance" = var,
    "Standard deviation" = sd,
    "Coefficient of Variation" = function(x, na.rm = FALSE) {
        sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm) },
    "Fano Factor" = function(x, na.rm = FALSE) {
        var(x, na.rm = na.rm) / mean(x, na.rm = na.rm) },
    "Mean" = mean,
    "Median" = median
)

# BACK-END =========================================================================================

server = function(input, output) {

    # Updated dataset selection based on cell type
    output$dataset_selection = renderUI({
        subt = md[md$cell_type == input$cell_type,]
        selectInput("dataset_label", "Choose a dataset",
            choices = sort(unique(paste0(subt$dataset, '/', subt$label))))
    })

    # Description of factor columns in md
    md_factor_descr = renderUI({
        lapply(colnames(md)[which(md_col_type == "factor")], FUN = function(name) {
            p(strong(md_col_label[[name]]), ":", md_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("md_factor_descr_", i)]] = md_factor_descr
    })

    # Description of real columns in md
    md_real_descr = renderUI({
        lapply(colnames(md)[which(md_col_type == "real")], FUN = function(name) {
            p(strong(md_col_label[[name]]), ":", md_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("md_real_descr_", i)]] = md_real_descr
    })

    # Description of factor columns in mda
    mda_factor_descr = renderUI({
        lapply(colnames(mda)[which(mda_col_type == "factor")], FUN = function(name) {
            p(strong(mda_col_label[[name]]), ":", mda_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("mda_factor_descr_", i)]] = mda_factor_descr
    })

    # Description of real columns in mda
    mda_real_descr = renderUI({
        lapply(colnames(mda)[which(mda_col_type == "real")], FUN = function(name) {
            p(strong(mda_col_label[[name]]), ":", mda_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("mda_real_descr_", i)]] = mda_real_descr
    })

    # Description of factor columns in nd
    nd_factor_descr = renderUI({
        lapply(colnames(nd)[which(nd_col_type == "factor")], FUN = function(name) {
            p(strong(nd_col_label[[name]]), ":", nd_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("nd_factor_descr_", i)]] = nd_factor_descr
    })

    # Description of real columns in nd
    nd_real_descr = renderUI({
        lapply(colnames(nd)[which(nd_col_type == "real")], FUN = function(name) {
            p(strong(nd_col_label[[name]]), ":", nd_col_descr[[name]])
        })
    })
    l = lapply(1:100, FUN = function(i) {
        output[[paste0("nd_real_descr_", i)]] = nd_real_descr
    })

    # Dataset --------------------------------------------------------------------------------------

    output$dataset_specs = renderUI({
        data = select_dot_data(input)
        out = list()
        out = list(out, h3("General"))

        # Cell type(s)
        cell_types = unique(data$cell_type)
        if( 1 != length(cell_types) ) {
            cell_types = lapply(unique(data$cell_type), FUN = function(x) span(paste0(x, ",")))
        }
        out = list(out, p(strong("Cell type"), ":", cell_types))

        # Dataset(s)
        datasets = unique(data$dataset)
        if( 1 != length(datasets) ) {
            datasets = lapply(unique(data$dataset), FUN = function(x) span(paste0(x, ",")))
        }
        out = list(out, p(strong("Dataset"), ":", datasets))

        # Number of cells
        out = list(out, h3("Cells"))
        ncell = paste0(data$dataset, "~", data$File, "~", data$cell_ID)[which(data$cell_ID != 0)]
        ncell = length(unique(ncell))
        out = list(out, p(strong("Number of cells"), ":", ncell))
        ncell = paste0(data$dataset, "~", data$File, "~", data$cell_ID)
        ncell = length(unique(ncell[which(data$cell_ID != 0 & data$G1 == 1)]))
        out = list(out, p(strong("Number of G1 cells"), ":", ncell))

        # Number of dots
        out = list(out, h3("Dots"))
        out = list(out, p(strong("Number of dots"), ":", nrow(data)))
        out = list(out, p(strong("Number of dots within cells"), ":",
            length(which(data$cell_ID != 0))))
        out = list(out, p(strong("Number of dots outside cells"), ":",
            length(which(data$cell_ID == 0))))
        out = list(out, p(strong("Number of allele couples"), ":",
            length(which(data$Allele == "Central"))))

        return(out)
    })

    output$dots_barPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_abs_xaxis, md_col_label)
        yax = input$dots_abs_yaxis

        # Factorize x-axis for color-scaling
        data[, xax] = factor(data[, xax])
        data[, xax] = order_data(data[, xax], xax)

        # Absolute number of dots
        if( yax == "All dots" ) {
            t = table(data[, xax])
            if ( input$dots_relative ) { t = t / table(data[, xax]) * 100 }
            t = data.frame(
                x = names(t),
                y = as.numeric(t)
            )
        }

        # Absolute number of dots in cells
        if( yax == "Dots in cells" ) {
            if( 0 != length(which(data$cell_ID != 0)) ) {
                t = table(data[data$cell_ID != 0, xax])
                if ( input$dots_relative ) {
                    t = t / table(data[, xax]) * 100
                }
                t = data.frame(
                    x = names(t),
                    y = as.numeric(t)
                )
            } else {
                t = data.frame(
                    x = unique(data[, xax]),
                    y = rep(0, length(unique(data[, xax])))
                )
            }
        }

        # Absolute number of dots outside cells
        if( yax == "Dots outside cells" ) {
            if( 0 != length(which(data$cell_ID == 0)) ) {
                t = table(data[data$cell_ID == 0, xax])
                if ( input$dots_relative ) {
                    t = t / table(data[, xax]) * 100
                }
                t = data.frame(
                    x = names(t),
                    y = as.numeric(t)
                )
            } else {
                t = data.frame(
                    x = unique(data[, xax]),
                    y = rep(0, length(unique(data[, xax])))
                )
            }
        }

        t$x = order_data(t$x, xax)
        p = ggplot(t, aes(x = x, y = y, fill = x))
        p = p + geom_bar(stat = "identity")
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

        if ( input$dots_relative ) {
            p = p + xlab(input$dots_abs_xaxis) + ylab(paste0(yax, " (%)"))
            p = p + ylim(0, 100)
        } else {
            p = p + xlab(input$dots_abs_xaxis) + ylab(yax)
        }

        p
    })

    output$nuclei_distPlot = renderPlotly({
        data = select_nuclei_data(input)

        # Identify axes
        xax = label2col(input$nuclei_dist_xaxis, nd_col_label)

        p = ggplot(data, aes_string(x = xax))
        if( input$nd_distr_line ) {
            p = p + geom_density(color = 1, fill = 'white')
        } else {
            p = p + geom_histogram(color = 1, fill = 'white', bins = input$nuclei_dist_nbins)
        }
        p = p + xlab(input$nuclei_dist_xaxis)

        # If only one dataset is selected, add G1 ranges
        if( 1 == length(unique(data$dataset)) ) {
            p = p + geom_vline(xintercept = min(data[data$G1 == 1, xax], na.rm = T),
                color = "red", linetype = "dashed")
            p = p + geom_vline(xintercept = max(data[data$G1 == 1, xax], na.rm = T),
                color = "red", linetype = "dashed")
        }

        p
    })

    output$nuclei_table = renderDataTable({
        data = select_nuclei_data(input)
        return(data[, c("cell_type", "dataset", "s", "n", "G1", "flat_size", "size", "surf",
            "sumI", "meanI", "shape")])
    })

    # Plot dots ------------------------------------------------------------------------------------

    output$dots_distPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify x axis
        xax = label2col(input$dots_dist_xaxis, md_col_label)

        p = ggplot(data, aes_string(x = xax))
        if( input$dots_dist_line ) {
            p = p + geom_density(color = 1, fill = 'white')
        } else {
            p = p + geom_histogram(color = 1, fill = 'white', bins = input$dots_dist_nbins)
        }
        p = p + xlab(input$dots_dist_xaxis)
    })

    output$dots_boxPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_box_xaxis, md_col_label)
        data[, xax] = factor(data[, xax])
        data[, xax] = order_data(data[, xax], xax)
        yax = label2col(input$dots_box_yaxis, md_col_label)

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax)) + geom_boxplot()
        p = p + xlab(input$dots_box_xaxis) + ylab(input$dots_box_yaxis)
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

        p
    })

    output$dots_sumPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_sum_xaxis, md_col_label)
        if ( "factor" == md_col_type[[xax]] ) {
            data[, xax] = factor(data[, xax])
            data[, xax] = order_data(data[, xax], xax)
        }
        yax = label2col(input$dots_sum_yaxis, md_col_label)
        stype = sum_types[[input$dots_sum_type]]

        # Remove missing x
        data = data[!is.na(data[, xax]),]

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax))
        p = p + stat_summary(fun.y = stype, geom = "point")
        p = p + xlab(input$dots_sum_xaxis)
        p = p + ylab(paste0(input$dots_sum_type, " of ", input$dots_sum_yaxis))
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

        if ( input$dots_sum_lm ) {

            if ( "factor" == md_col_type[[xax]] ) {
                sdata = data.frame(
                    x = levels(data[, xax]),
                    y = as.numeric(by(data[, yax], data[, xax], FUN = stype, na.rm = T))
                )
                sdata = sdata[!is.na(sdata$y),]
                regr = summary(lm(sdata$y ~ I(1:nrow(sdata))))$coefficients
            } else {
                sdata = do.call(rbind, by(data, data[, xax],
                    FUN = function(subt) {
                        data.frame(
                            x = subt[1, xax],
                            y = stype(subt[, yax], na.rm = T)
                        )
                    }
                ))
                sdata = sdata[!is.na(sdata$y),]
                regr = summary(lm(sdata$y ~ sdata$x))$coefficients
            }
            p = p + geom_abline(intercept = regr[1], slope = regr[2], linetype = 2)
        }

        p
    })

    output$dots_summary_lm = renderPrint({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_sum_xaxis, md_col_label)
        if ( "factor" == md_col_type[[xax]] ) {
            data[, xax] = factor(data[, xax])
            data[, xax] = order_data(data[, xax], xax)
        }
        yax = label2col(input$dots_sum_yaxis, md_col_label)
        stype = sum_types[[input$dots_sum_type]]

        # Remove missing x
        data = data[!is.na(data[, xax]),]

        if ( "factor" == md_col_type[[xax]] ) {
            sdata = data.frame(
                x = levels(data[, xax]),
                y = as.numeric(by(data[, yax], data[, xax], FUN = stype, na.rm = T))
            )
            sdata = sdata[!is.na(sdata$y),]
            regr = summary(lm(sdata$y ~ I(1:nrow(sdata))))
        } else {
            sdata = do.call(rbind, by(data, data[, xax],
                FUN = function(subt) {
                    data.frame(
                        x = subt[1, xax],
                        y = stype(subt[, yax], na.rm = T)
                    )
                }
            ))
            sdata = sdata[!is.na(sdata$y),]
            regr = summary(lm(sdata$y ~ sdata$x))
        }

        return(regr)
    })

    output$dots_scatterPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_scat_xaxis, md_col_label)
        yax = label2col(input$dots_scat_yaxis, md_col_label)
        col = label2col(input$dots_scat_color, md_col_label)
        data[, col] = factor(data[, col])
        data[, col] = order_data(data[, col], col)

        p = ggplot(data, aes_string(x = xax, y = yax, color = col)) + geom_point()
        p = p + xlab(input$dots_scat_xaxis) + ylab(input$dots_scat_yaxis)

        if ( input$dots_scat_lm ) {
            regr = summary(lm(data[, yax] ~ data[, xax]))$coefficients
            p = p + geom_abline(intercept = regr[1], slope = regr[2], linetype = 2)
        }

        p
    })

    output$dots_scatter_lm = renderPrint({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_scat_xaxis, md_col_label)
        yax = label2col(input$dots_scat_yaxis, md_col_label)
        regr = summary(lm(data[, yax] ~ data[, xax]))

        return(regr)
    })

    output$dots_hexPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify axes
        xax = label2col(input$dots_hex_xaxis, md_col_label)
        yax = label2col(input$dots_hex_yaxis, md_col_label)

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax)) + geom_hex()
        p = p + xlab(input$dots_hex_xaxis) + ylab(input$dots_hex_yaxis)
    })

    output$single_dot_table = renderDataTable({
        data = select_dot_data(input)
        return(data[, c("cell_type", "dataset", "File", "label", "chr", "set", "Channel", "cell_ID",
            "G1", "x", "y", "z", "lamin_dist", "lamin_dist_norm", "centr_dist", "centr_dist_norm",
            "Allele")])
    })

    # Plot allele ----------------------------------------------------------------------------------

    output$allele_distPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify x axis
        xax = label2col(input$allele_dist_xaxis, mda_col_label)

        p = ggplot(data, aes_string(x = xax))
        if( input$allele_dist_line ) {
            p = p + geom_density(color = 1, fill = 'white')
        } else {
            p = p + geom_histogram(color = 1, fill = 'white', bins = input$allele_dist_nbins)
        }
        p = p + xlab(input$allele_dist_xaxis)

        p
    })

    output$allele_boxPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_box_xaxis, mda_col_label)
        data[, xax] = factor(data[, xax])
        data[, xax] = order_data(data[, xax], xax)
        yax = label2col(input$allele_box_yaxis, mda_col_label)

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax)) + geom_boxplot()
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        p = p + xlab(input$allele_box_xaxis) + ylab(input$allele_box_yaxis)

        p
    })

    output$allele_sumPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_sum_xaxis, mda_col_label)
        if ( "factor" == mda_col_type[[xax]] ) {
            data[, xax] = factor(data[, xax])
            data[, xax] = order_data(data[, xax], xax)
        }
        yax = label2col(input$allele_sum_yaxis, mda_col_label)
        stype = sum_types[[input$allele_sum_type]]

        # Remove missing x
        data = data[!is.na(data[, xax]),]

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax))
        p = p + stat_summary(fun.y = stype, geom = "point")
        p = p + xlab(input$allele_box_xaxis)
        p = p + ylab(paste0(input$allele_sum_type, " of ", input$allele_box_yaxis))
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

        if ( input$allele_sum_lm ) {
            if ( "factor" == mda_col_type[[xax]] ) {
                sdata = data.frame(
                    x = levels(data[, xax]),
                    y = as.numeric(by(data[, yax], data[, xax], FUN = stype, na.rm = T))
                )
                sdata = sdata[!is.na(sdata$y),]
                regr = summary(lm(sdata$y ~ I(1:nrow(sdata))))$coefficients
            } else {
                sdata = do.call(rbind, by(data, data[, xax],
                    FUN = function(subt) {
                        data.frame(
                            x = subt[1, xax],
                            y = stype(subt[, yax], na.rm = T)
                        )
                    }
                ))
                sdata = sdata[!is.na(sdata$y),]
                regr = summary(lm(sdata$y ~ sdata$x))$coefficients
            }
            p = p + geom_abline(intercept = regr[1], slope = regr[2], linetype = 2)
        }

        p
    })

    output$allele_summary_lm = renderPrint({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_sum_xaxis, mda_col_label)
        if ( "factor" == mda_col_type[[xax]] ) {
            data[, xax] = factor(data[, xax])
            data[, xax] = order_data(data[, xax], xax)
        }
        yax = label2col(input$allele_sum_yaxis, mda_col_label)
        stype = sum_types[[input$allele_sum_type]]

        # Remove missing x
        data = data[!is.na(data[, xax]),]

        if ( "factor" == mda_col_type[[xax]] ) {
            sdata = data.frame(
                x = levels(data[, xax]),
                y = as.numeric(by(data[, yax], data[, xax], FUN = stype, na.rm = T))
            )
            sdata = sdata[!is.na(sdata$y),]
            regr = summary(lm(sdata$y ~ I(1:nrow(sdata))))
        } else {
            sdata = do.call(rbind, by(data, data[, xax],
                FUN = function(subt) {
                    data.frame(
                        x = subt[1, xax],
                        y = stype(subt[, yax], na.rm = T)
                    )
                }
            ))
            sdata = sdata[!is.na(sdata$y),]
            regr = summary(lm(sdata$y ~ sdata$x))
        }

        return(regr)
    })

    output$allele_scatterPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_scat_xaxis, mda_col_label)
        yax = label2col(input$allele_scat_yaxis, mda_col_label)
        col = label2col(input$allele_scat_color, mda_col_label)
        data[, col] = factor(data[, col])
        data[, col] = order_data(data[, col], col)

        p = ggplot(data, aes_string(x = xax, y = yax, color = col)) + geom_point()
        p = p + xlab(input$allele_scat_xaxis) + ylab(input$allele_scat_yaxis)

        if ( input$allele_scat_lm ) {
            regr = summary(lm(data[, yax] ~ data[, xax]))$coefficients
            p = p + geom_abline(intercept = regr[1], slope = regr[2], linetype = 2)
        }

        p
    })

    output$allele_scatter_lm = renderPrint({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_scat_xaxis, mda_col_label)
        yax = label2col(input$allele_scat_yaxis, mda_col_label)
        regr = summary(lm(data[, yax] ~ data[, xax]))

        return(regr)
    })

    output$allele_hexPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify axes
        xax = label2col(input$allele_hex_xaxis, mda_col_label)
        yax = label2col(input$allele_hex_yaxis, mda_col_label)

        p = ggplot(data, aes_string(x = xax, y = yax, color = xax)) + geom_hex()
        p = p + xlab(input$allele_hex_xaxis) + ylab(input$allele_hex_yaxis)

        p
    })

    output$allele_table = renderDataTable({
        return(mda[, c("cell_type", "dataset", "File", "label", "chr", "set", "Channel", "cell_ID",
            "G1", "d_3d", "d_lamin", "d_lamin_norm", "d_centr", "d_centr_norm", "angle")])
    })

    # Rankings -------------------------------------------------------------------------------------

    output$rankPlot = renderPlotly({
        data = select_dot_data(input)

        # Identify x axis
        xax = label2col(input$rank_xaxis, md_col_label)
        data[, xax] = factor(data[, xax])
        data[, xax] = order_data(data[, xax], xax)
        yax = label2col(input$rank_yaxis, md_col_label)

        sdata = do.call(rbind, by(data, data[, xax],
            FUN = function(subt) {
                data.frame(
                    x = subt[1, xax],
                    y = sum_types[[input$rank_stype]](subt[, yax], na.rm = T)
                )
            }
        ))
        sdata$x = factor(sdata$x, levels = sdata$x[order(sdata$y)])

        p = ggplot(sdata, aes(x, y, fill = x))
        p = p + geom_bar(stat = "identity")
        p = p + xlab(input$rank_xaxis) + ylab(paste0(input$rank_stype, " of ", input$rank_yaxis))
        p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

        p
    })

    output$rank_table = renderPrint({
       data = select_dot_data(input)

        # Identify x axis
        xax = label2col(input$rank_xaxis, md_col_label)
        data[, xax] = factor(data[, xax])
        data[, xax] = order_data(data[, xax], xax)
        yax = label2col(input$rank_yaxis, md_col_label)

        sdata = do.call(rbind, by(data, data[, xax],
            FUN = function(subt) {
                data.frame(
                    x = subt[1, xax],
                    y = sum_types[[input$rank_stype]](subt[, yax], na.rm = T)
                )
            }
        ))
        sdata$x = factor(sdata$x, levels = sdata$x[order(sdata$y)])
        rownames(sdata) = NULL
        colnames(sdata) = c("chr", input$rank_stype)

        return(sdata[order(sdata[,2]),])
    })

    # E/C dot sets plots ---------------------------------------------------------------------------
    
    output$ecset_distPlot = renderPlotly({
        data = select_allele_data(input)

        # Identify x axis
        xax = label2col(input$ecset_dist_xaxis, mda_col_label)

        p = ggplot(data, aes_string(x = xax, fill = "set", colour = "set"))
        p = p + geom_density(fill = NA)
        p = p + xlab(input$ecset_dist_xaxis)

        p
    })

}

# FRONT-END ========================================================================================

ui <- fluidPage(

    # App title
    titlePanel("GPSeq/iFISH project"),

    helpText(paste0("Dataset generated on ", dataset_date, ".")),

    # Sidebar panel for inputs ---------------------------------------------------------------------
    sidebarPanel(
        conditionalPanel("!input.all_cell_types",
            selectInput("cell_type", label = "Choose a cell type",
                choices = unique(md$cell_type), selected = "IMR90")
        ),
        checkboxInput("all_cell_types", "include all cell types", FALSE),
        conditionalPanel("!input.all_datasets", uiOutput("dataset_selection")),
        checkboxInput("all_datasets", "include all datasets", FALSE),
        selectInput("cell_phase", label = "Choose a cell sub-population",
            choices = c(
                "All cells",
                "Only cells with 1 allele",
                "Only cells with 2 alleles"
            ), selected = "IMR90"),
        checkboxInput("g1only", "use only G1 cells.", FALSE)
    ),

    # Main panel
    mainPanel(tabsetPanel(

        # General view -----------------------------------------------------------------------------
        tabPanel("General", br(),
            tabsetPanel(
                tabPanel("Home",

                    h3("Introduction"),

                    p("
                        This interface allows to navigate GPSeq/iFISH data.
                        The purpose is to navigate the data, and identify possible trends/patterns
                        or new ways of plotting the data to highlight such trends/patterns.
                    "),

                    p("Five ", strong("main tabs"), " are available:"),
                    tags$ul(
                        tags$li(
                            strong("General"),
                            span(": general remarks, instructions and help page.")
                        ),
                        tags$li(
                            strong("Dataset"),
                            span(": dataset(s) features.")
                        ),
                        tags$li(
                            strong("Single-dot"), span(": single-dot data.")
                        ),
                        tags$li(
                            strong("Allele pairs"),
                            span(": data on dots couple per channel, considered as alleles.")
                        ),
                        tags$li(
                            strong("Specific"),
                            span(": specific plots requested by the users,
                                not achievable otherwise through the interface.")
                        )
                    ),

                    p("Each tab has a number of ", strong("sub-tabs"), ", which correspond to
                        different types of plot. Each plot is plotted using", strong("plotly"),
                        "which allows to interact with the plots (zoom, pan, see single-point
                        values). At the same time, the intearface contains a number of interactive
                        fields to select different variables as X/Y axes, grouping etc..."),

                    p("More information are available in each sub-tab and in the help page.")

                ),
                tabPanel("Help",
                    br(),
                    p("Lorem ipsum dolor sit amet, consectetur adipisicing elit. Officia dolorem
                        ducimus veritatis dolores assumenda, iusto qui porro nihil vitae expedita
                        aperiam eos magnam totam officiis, accusamus asperiores esse a ipsam.")
                ),
                tabPanel("To-do",
                    br(),
                    tags$ul(
                        tags$li(tags$s("Add option for linear model fitting to scatter plots.")),
                        tags$li(tags$s("Add option to plot the variability of
                            a continuous variable against a categorical one.")),
                        tags$li("Add plots to specific tab.")
                    )
                )
            )
        ),

        # Nuclear data -----------------------------------------------------------------------------
        tabPanel("Dataset",
            br(),

            tabsetPanel(
                tabPanel("Specs",
                    br(), uiOutput(outputId = "dataset_specs")
                ),
                tabPanel("Dots",
                    plotlyOutput(outputId = "dots_barPlot"), br(),
                    fluidRow(column(6,
                        checkboxInput("dots_relative", "show percentage of dots.", FALSE)
                    )),
                    fluidRow(
                        column(6,
                            selectInput("dots_abs_yaxis", label = "Y-axis",
                                choices = c("All dots", "Dots in cells", "Dots outside cells")),

                                p(
                                    strong("All dots"),
                                    ": consider ", tags$u("all"), " dots."
                                ),
                                p(
                                    strong("Dots in cells"),
                                    ": consider only dots ", tags$u("within"), " cells."
                                ),
                                p(
                                    strong("Dots outside cells"),
                                    ": consider only dots ", tags$u("outside"), " cells."
                                )

                        ),
                        column(6,
                            selectInput("dots_abs_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "factor")],
                                    FUN = function(x) { md_col_label[[x]] }),
                                selected = "Channel"),
                            uiOutput("md_factor_descr_1")
                        )
                    )
                ),
                tabPanel("Nuclei distribution", br(),
                    h4(strong("NOTE:"), "cell sub-population filters do not apply to this plot.",
                        "Only the 'use only G1 cells' option can be used in this context."), br(),
                    plotlyOutput(outputId = "nuclei_distPlot"), br(),
                    fluidRow(
                        column(6,
                            checkboxInput("nd_distr_line", "replace histogram with density.", FALSE)
                        ),
                        column(6,
                            sliderInput("nuclei_dist_nbins", "Number of bins:",
                                min = 1, max = 300, value = 30)
                        )
                    ),
                    fluidRow(
                        column(6,
                            selectInput("nuclei_dist_xaxis", label = "X-axis",
                                choices = lapply(colnames(nd)[which(nd_col_type == "real")],
                                    FUN = function(x) { nd_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("nd_real_descr_1")
                        )
                    )
                ),
                tabPanel("Raw table", br(),
                    dataTableOutput(outputId = "nuclei_table")
                )
            )
        ),

        # Single dot data --------------------------------------------------------------------------
        tabPanel("Single-dot",
            br(), p("Here you can find the data on single FISH dots."),

            tabsetPanel(
                tabPanel("Distribution",
                    plotlyOutput(outputId = "dots_distPlot"), br(),
                    fluidRow(
                        column(6,
                            checkboxInput("dots_dist_line",
                                "replace histogram with density.", FALSE)
                        ),
                        column(6,
                            sliderInput("dots_dist_nbins", "Number of bins:",
                                min = 1, max = 100, value = 30)
                        )
                    ),
                    fluidRow(
                        column(6,
                            selectInput("dots_dist_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("md_real_descr_1")
                        )
                    )
                ),
                tabPanel("Boxplot",
                    plotlyOutput(outputId = "dots_boxPlot"), br(),
                    fluidRow(
                        column(6,
                            selectInput("dots_box_yaxis", label = "Y-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("md_real_descr_2")
                        ),
                        column(6,
                            selectInput("dots_box_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "factor")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Allele label (by GPSeq)"),
                            uiOutput("md_factor_descr_2")
                        )
                    )
                ),
                tabPanel("Summary",
                    plotlyOutput(outputId = "dots_sumPlot"), br(),
                    fluidRow(
                        column(12,
                            checkboxInput("dots_sum_lm", "perform linear regression.", FALSE),
                            conditionalPanel("input.dots_sum_lm",
                                verbatimTextOutput("dots_summary_lm")
                            )
                        )
                    ),
                    fluidRow(
                        column(6,
                            selectInput("dots_sum_xaxis", label = "X-axis",
                                choices = c(
                                    lapply(colnames(md)[which(md_col_type == "factor")],
                                        FUN = function(x) { md_col_label[[x]] }),
                                    md_col_label[["chr_size"]]),
                                selected = "Allele label (by GPSeq)"),
                            uiOutput("md_factor_descr_4")
                        ),
                        column(6,
                            selectInput("dots_sum_type", label = "Summary measure (Y-axis)",
                                choices = names(sum_types), selected = "Variance"),
                            selectInput("dots_sum_yaxis", label = "Y-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("md_real_descr_4")
                        )
                    )
                ),
                tabPanel("Scatter",
                    plotlyOutput(outputId = "dots_scatterPlot"), br(),
                    fluidRow(
                        column(12,
                            checkboxInput("dots_scat_lm", "perform linear regression.", FALSE),
                            conditionalPanel("input.dots_scat_lm",
                                verbatimTextOutput("dots_scatter_lm")
                            )
                        )
                    ),
                    fluidRow(
                        column(6,
                            selectInput("dots_scat_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            selectInput("dots_scat_yaxis", label = "Y-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Relative distance from lamina [a.u.]"),
                            uiOutput("md_real_descr_3")
                        ),
                        column(6,
                            selectInput("dots_scat_color", label = "Color",
                                choices = lapply(colnames(md)[which(md_col_type == "factor")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Allele label (by GPSeq)"),
                            uiOutput("md_factor_descr_3")
                        )
                    )
                ),
                tabPanel("Hexplot",
                    plotlyOutput(outputId = "dots_hexPlot"), br(),
                    fluidRow(
                        column(6,
                            selectInput("dots_hex_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("md_real_descr_5")
                        ),
                        column(6,
                            selectInput("dots_hex_yaxis", label = "Y-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Relative distance from lamina [a.u.]"),
                            uiOutput("md_real_descr_6")
                        )
                    )
                ),
                tabPanel("Raw table", br(),
                    dataTableOutput(outputId = "single_dot_table")
                )
            )
        ),

        # Allele data ------------------------------------------------------------------------------
        tabPanel("Allele pairs",
            br(), p("Here you can find the data on allele couples; thus,",
                strong("only cells with 2 alleles are considered here.")),

            conditionalPanel("input.cell_phase != 'Only cells with 1 allele'",
                tabsetPanel(
                    tabPanel("Distribution",
                        plotlyOutput(outputId = "allele_distPlot"), br(),
                        fluidRow(
                            column(6,
                                checkboxInput("allele_dist_line",
                                    "replace histogram with density.", FALSE)
                            ),
                            column(6,
                                sliderInput("allele_dist_nbins", "Number of bins:",
                                    min = 1, max = 100, value = 30)
                            )
                        ),
                        fluidRow(
                            column(6,
                                selectInput("allele_dist_xaxis", label = "X-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Absolute 3D distance [nm]"),
                                uiOutput("mda_real_descr_1")
                            )
                        )
                    ),
                    tabPanel("Boxplot",
                        plotlyOutput(outputId = "allele_boxPlot"), br(),
                        fluidRow(
                            column(6,
                                selectInput("allele_box_yaxis", label = "Y-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Absolute 3D distance [nm]"),
                                uiOutput("mda_real_descr_2")
                            ),
                            column(6,
                                selectInput("allele_box_xaxis", label = "X-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "factor")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Channel"),
                                uiOutput("mda_factor_descr_1")
                            )
                        )
                    ),
                    tabPanel("Summary",
                        plotlyOutput(outputId = "allele_sumPlot"), br(),
                        fluidRow(
                            column(12,
                                checkboxInput("allele_sum_lm", "perform linear regression.", FALSE),
                                conditionalPanel("input.allele_sum_lm",
                                    verbatimTextOutput("allele_summary_lm")
                                )
                            )
                        ),
                        fluidRow(
                            column(6,
                                selectInput("allele_sum_xaxis", label = "X-axis",
                                    choices = c(
                                        lapply(colnames(mda)[which(mda_col_type == "factor")],
                                            FUN = function(x) { mda_col_label[[x]] }),
                                        "Chromosome length [bp]"
                                    ), selected = "Channel"),
                                uiOutput("mda_factor_descr_4")
                            ),
                            column(6,
                                selectInput("allele_sum_type", label = "Summary measure (Y-axis)",
                                    choices = names(sum_types), selected = "Variance"),
                                selectInput("allele_sum_yaxis", label = "Y-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Absolute 3D distance [nm]"),
                                uiOutput("mda_real_descr_4")
                            )
                        )
                    ),
                    tabPanel("Scatter",
                        plotlyOutput(outputId = "allele_scatterPlot"), br(),
                        fluidRow(
                            column(12,
                                checkboxInput("allele_scat_lm",
                                    "perform linear regression.", FALSE),
                                conditionalPanel("input.allele_scat_lm",
                                    verbatimTextOutput("allele_scatter_lm")
                                )
                            )
                        ),
                        fluidRow(
                            column(6,
                                selectInput("allele_scat_xaxis", label = "X-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Absolute 3D distance [nm]"),
                                selectInput("allele_scat_yaxis", label = "Y-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Distance between normalized lamina layers [a.u.]"),
                                uiOutput("mda_real_descr_3")
                            ),
                            column(6,
                                selectInput("allele_scat_color", label = "Color",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "factor")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Channel"),
                                uiOutput("mda_factor_descr_2")
                            )
                        )
                    ),
                    tabPanel("Hexplot",
                        plotlyOutput(outputId = "allele_hexPlot"), br(),
                        fluidRow(
                            column(6,
                                selectInput("allele_hex_xaxis", label = "X-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Absolute 3D distance [nm]"),
                                uiOutput("mda_real_descr_5")
                            ),
                            column(6,
                                selectInput("allele_hex_yaxis", label = "Y-axis",
                                    choices = lapply(colnames(mda)[which(mda_col_type == "real")],
                                    FUN = function(x) { mda_col_label[[x]] }),
                                    selected = "Distance between normalized lamina layers [a.u.]"),
                                uiOutput("mda_real_descr_6")
                            )
                        )
                    ),
                    tabPanel("Raw table", br(),
                        dataTableOutput(outputId = "allele_table")
                    )
                )
            ),
            conditionalPanel("input.cell_phase == 'Only cells with 1 allele'",
                br(),
                p("Allele couple data are not available when selecting only cells with one allele.")
            )
        ),

        # Specific plots ---------------------------------------------------------------------------
        tabPanel("Specific",
            br(), p("Here you can specific user-requested plots and features."),

            tabsetPanel(
                tabPanel("Rankings",
                    br(),
                    p("Here to plot chromosome rankings based on centrality."),
                    plotlyOutput(outputId = "rankPlot"), br(),
                    conditionalPanel("input.rank_table",
                        verbatimTextOutput("rank_table")
                    ),
                    fluidRow(
                        column(6,
                            checkboxInput("rank_table", "show rank table.", FALSE)
                        ),
                        column(6,
                            selectInput("rank_stype", label = "Summary measure",
                                choices = names(sum_types), selected = "Mean")
                        )
                    ),
                    fluidRow(
                        column(6,
                            selectInput("rank_xaxis", label = "X-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "factor")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Chromosome"),
                            uiOutput("md_factor_descr_8")
                        ),
                        column(6,
                            selectInput("rank_yaxis", label = "Y-axis",
                                choices = lapply(colnames(md)[which(md_col_type == "real")],
                                FUN = function(x) { md_col_label[[x]] }),
                                selected = "Absolute distance from lamina [nm]"),
                            uiOutput("md_real_descr_8")
                        )
                    )
                ),
                tabPanel("E/C-sets",
                    br(),
                    p("Here characterization of e/c probe sets.",
                        "'e' probes are 'terminal', close to telomeres.",
                        "'c' probes are 'central', further from telomeres."),
                    
                    br(), p("Here you can find the data on allele couples; thus,",
                        strong("only cells with 2 alleles are considered here."),
                        "Moreover, here alleles are considered at the single-probe level."),

                    conditionalPanel("input.cell_phase != 'Only cells with 1 allele'",
                        tabsetPanel(
                            tabPanel("Overall distribution",
                                plotlyOutput(outputId = "ecset_distPlot"), br(),
                                fluidRow(
                                ),
                                fluidRow(
                                    column(6,
                                        selectInput("ecset_dist_xaxis", label = "X-axis",
                                            choices = lapply(colnames(mda)
                                                [which(mda_col_type == "real")],
                                            FUN = function(x) { mda_col_label[[x]] }),
                                            selected = "Absolute 3D distance [nm]"),
                                        uiOutput("mda_real_descr_9")
                                    )
                                )
                            )
                        )
                    )
                ),
                tabPanel("Linear/Spatial",
                    br(),
                    p("Here to plot linear genomic distance against 3D spatial distance."),
                    h1("@TODO")
                )
            )
        )

    , selected = "Specific"))
)

# RUN ==============================================================================================

shinyApp(ui = ui, server = server)

# END ----------------------------------------------------------------------------------------------

####################################################################################################
