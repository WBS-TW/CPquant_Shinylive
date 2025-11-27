# This version was copied 20250910
# Uses Shinylive to deploy static web app
# see: https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html

library(ggplot2)
library(shiny)
library(htmlwidgets)
library(readxl)
library(stringr)
library(nnls)
library(dplyr)
library(tibble)
library(tidyr)
library(DT)
library(plotly)
library(purrr)
library(markdown)
library(openxlsx)
library(S7)


# Various utilities and helper functions for CPquant


#########################################################################################################
#------------------------------------------ CPquant_utils.R ------------------------------------------#
#########################################################################################################

### Function to perform deconvolution on a single data frame ###
perform_deconvolution <- function(df, combined_standard, CPs_standards_sum_RF) {

    df_matrix <- as.matrix(df)

    print(paste("df_matrix dimensions:", dim(df_matrix)))
    print(paste("combined_standard dimensions:", dim(combined_standard)))

    if (nrow(combined_standard) != nrow(df_matrix)) {
        stop("Dimensions of combined_standard and df are incompatible.")
    }

    # Reshape df_matrix if it has only one column or extract the first column if it has multiple
    if (ncol(df_matrix) == 1) {
        df_vector <- as.vector(df_matrix)
    } else {
        df_vector <- as.vector(df_matrix["Relative_Area"])  # Extract the first column for nnls
    }

    # Check for NA/NaN/Inf values in df_vector and combined_standard
    if (any(is.na(df_vector)) || any(is.nan(df_vector)) || any(is.infinite(df_vector))) {
        stop("df_vector contains NA/NaN/Inf values.")
    }

    if (any(is.na(combined_standard)) || any(is.nan(combined_standard)) || any(is.infinite(combined_standard))) {
        stop("combined_standard contains NA/NaN/Inf values.")
    }

    # Perform nnls
    # combined_standard is response factor and df_vector is Relative_Area
    deconv <- nnls::nnls(combined_standard, df_vector)

    # Extract deconvolution results
    deconv_coef <- deconv$x


    #Normalize the coefficients so they sum to 1
    if (sum(deconv_coef) > 0) {
        deconv_coef <- deconv_coef / sum(deconv_coef)
    }else(deconv_coef <- 0)

    # Calculate deconvolved values (which is the response factor) using matrix multiplication
    deconv_resolved <- combined_standard %*% deconv_coef

    # Calculate the sum of RF*frac
    sum_deconv_RF <- as.matrix(CPs_standards_sum_RF) %*% deconv_coef





    # Calculate the goodness of fit by coefficient of determination (R2)
    # Calculate the total sum of squares (SST)
    sst <- sum((df_vector - mean(df_vector))^2)
    # Calculate the residual sum of squares (SSR)
    ssr <- sum((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate R-squared coefficient of determination
    deconv_rsquared <- 1 - (ssr / sst)

    # Calculate Mean Squared Error (MSE)
    mse <- mean((df_vector - deconv_resolved/sum(deconv_resolved))^2)

    # Calculate Root Mean Squared Error (RMSE)
    rmse <- sqrt(mse)

    #Chiq-square test: ensure that values are positive for chi-square test
    # if (any(deconv_resolved < 0) || any(df_vector < 0)) {
    #     warning("Non-positive values found, skipping chi-square test")
    #     chisq_result <- NULL
    # } else {
    #     #adding a very small constant to avoid 0 values
    #     observed_corr <- df_vector + 1E-12
    #     predicted_corr <- deconv_resolved + 1E-12
    #     chisq_result <- chisq.test(x= observed_corr, p = predicted_corr/sum(predicted_corr), rescale.p = TRUE)
    # }


    # Kolmogorov-Smirnov Test
    #ks_result <- ks.test(deconv_resolved, df_vector)



    #combine results for output
    combined_standard_names <- colnames(combined_standard)

    names(deconv_coef) <- combined_standard_names


    return(list(
        sum_deconv_RF = sum_deconv_RF,
        deconv_coef = deconv_coef,
        deconv_resolved = deconv_resolved,
        deconv_rsquared = deconv_rsquared
        #chisq_result = chisq_result
    ))
}




################################################################################
##--------------------------- CPquant UI Components --------------------------##
################################################################################


defineVariablesUI <- function(Skyline_output){
    ###START: Define UI components

    # Create the UI components
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "removeSamples", #select if some samples will be removed from quantification
                label = 'Remove samples from quantification?',
                choices = unique(Skyline_output$Replicate_Name),
                selected = NULL,
                multiple = TRUE
            )
        ),

        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
        shiny::column(
            6,
            shiny::sliderInput(
                inputId = "removeRsquared", #keep only Molecule from standard calibration curves above this rsquared, zero means keep everything
                label = 'Keep the the calibration curves above this rsquared (0 means keep everything)',
                min = 0,
                max = 1,
                value = 0.70,
                step = 0.05
            )
            # shiny::column(
            #     6,
            #     shiny::checkboxInput(
            #           inputId = "zerointercept", #force y-intercept through zero
            #           label = "Set intercept to zero",
            #           value = FALSE
            #     )
        )
    )
}
### END FUNCTION


defineCorrectionUI <- function(Skyline_output){
    ###START: Define UI components

    # Create the UI components
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "chooseRS", #select which will be the RS
                label = 'Choose RS for correction',
                choices = unique(Skyline_output$Molecule[Skyline_output$Molecule_List == "RS"]),
                selected = NULL,
                multiple = FALSE
            )
        ))
}
### END FUNCTION


defineCalcrecoveryUI <- function(Skyline_output){
    ###START: Define UI components

    # Create the UI components
    shiny::fluidRow(
        shiny::column(
            6,
            shiny::selectInput(
                inputId = "chooseRS2", #select which will be the RS
                label = 'Choose RS for recovery',
                choices = unique(Skyline_output$Molecule[Skyline_output$Molecule_List == "RS"]),
                selected = NULL,
                multiple = FALSE
            )
        ))
}
### END FUNCTION


################################################################################
##------------------------------ CPquant plots -------------------------------##
################################################################################
#############################################################################
#############################################################################
plot_skyline_output <- function(Skyline_output){

    Skyline_output |>
        dplyr::filter(Isotope_Label_Type == "Quan") |>
        dplyr::mutate(OrderedMolecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)]))) |>  # Create a composite ordering factor
        plotly::plot_ly(
            x = ~ OrderedMolecule,
            y = ~ Area,
            color = ~ Sample_Type,
            type = "box",
            text = ~paste(
                "Homologue: ", PCA,
                "<br>Sample: ", Replicate_Name,
                "<br>Area:", round(Area, 2)
            ),
            hoverinfo = "text"
        ) |>
        plotly::layout(xaxis = list(title = 'Molecule'),
                       yaxis = list(title = 'Area'))
}

#############################################################################

plot_calibration_curves <- function(CPs_standards, quantUnit) {
    # Unnest the data first
    CPs_standards_unnested <- CPs_standards |>
        dplyr::filter(RF > 0) |>
        tidyr::unnest(data) |>
        dplyr::mutate(Molecule = factor(Molecule,
                                        levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique Quantification Groups
    groups <- unique(CPs_standards_unnested$Quantification_Group)

    # Create a list to store individual plots
    plot_list <- list()

    # Create individual plots for each group
    for(i in seq_along(groups)) {
        group_data <- CPs_standards_unnested |>
            dplyr::filter(Quantification_Group == groups[i])

        plot_list[[i]] <- plotly::plot_ly() |>
            plotly::add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~Area,
                color = ~PCA,
                type = 'scatter',
                mode = 'markers',
                legendgroup = ~Molecule_List,
                legendgrouptitle = list(text = ~Molecule_List),
                showlegend = FALSE,
                name = ~paste(PCA, "(points)"),
                text = ~paste(
                    "Molecule:", Molecule,
                    "<br>Area:", round(Area, 2),
                    "<br>Concentration:", round(Analyte_Concentration, 2),
                    "<br>R2:", round(cal_rsquared, 3)
                ),
                hoverinfo = 'text'
            ) |>
            plotly::add_trace(
                data = group_data,
                x = ~Analyte_Concentration,
                y = ~RF * Analyte_Concentration + intercept,
                color = ~PCA,
                type = 'scatter',
                mode = 'lines',
                legendgroup = ~Molecule_List,
                name = ~PCA,
                hoverinfo = 'none'
            )
    }

    # Calculate layout
    subplot_cols <- min(length(groups), 2)  # Maximum 2 columns
    subplot_rows <- ceiling(length(groups) / subplot_cols)

    # Create annotations for titles
    annotations <- list()
    for(i in seq_along(groups)) {
        row <- ceiling(i/subplot_cols)
        col <- if(i %% subplot_cols == 0) subplot_cols else i %% subplot_cols

        annotations[[i]] <- list(
            text = groups[i],
            font = list(size = 14),
            xref = "paper",
            yref = "paper",
            x = (col - 0.5)/subplot_cols,
            y = 1 - (row - 1)/subplot_rows,
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE
        )
    }

    # Combine plots using subplot
    final_plot <- plotly::subplot(
        plot_list,
        nrows = subplot_rows,
        shareX = FALSE,
        #shareY = TRUE,
        shareY = FALSE,
        margin = 0.1
    ) |>
        plotly::layout(
            height = 400 * subplot_rows,
            showlegend = TRUE,
            annotations = annotations,
            margin = list(t = 50, b = 50, l = 50, r = 50),
            legend = list(
                groupclick = "togglegroup",  # Changed from "toggleitem" to "togglegroup"
                tracegroupgap = 10,
                itemsizing = "constant"
            ),
            xaxis = list(title = paste0("Analyte Concentration/Amount (", quantUnit, ")"))  # Use the reactive value for x-axis label
        )

    return(final_plot)
}

#############################################################################

plot_quanqualratio <- function(Skyline_output_filt) {

    Skyline_output_filt |>
        dplyr::group_by(Replicate_Name, Molecule) |>
        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
        tidyr::fill(Quan_Area, .direction = "downup") |>
        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
        tidyr::fill(QuanMZ, .direction = "downup") |>
        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
        tidyr::replace_na(list(QuanQualRatio = 0)) |>
        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>
        dplyr::ungroup() |>
        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio) |>
        plotly::plot_ly(x = ~Replicate_Name, y = ~QuanQualRatio, type = 'violin', color = ~Sample_Type,
                        text = ~paste("Sample: ", Replicate_Name,
                                      "<br>Molecule List: ", Molecule_List,
                                      "<br>Molecule: ", Molecule,
                                      "<br>Quan/Qual MZ: ", QuanQualMZ,
                                      "<br>Ratio: ", round(QuanQualRatio, 2)),
                        hoverinfo = "text") |>
        plotly::layout(title = 'Quan-to-Qual Ratio',
                       xaxis = list(title = 'Replicate Name'),
                       yaxis = list(title = 'Quan-to-Qual Ratio'))
}

##############################################################################

plot_meas_vs_theor_ratio <- function(Skyline_output_filt) {

    Skyline_output_filt |>
        dplyr::group_by(Replicate_Name, Molecule) |>
        dplyr::mutate(Quan_Area = ifelse(Isotope_Label_Type == "Quan", Area, NA)) |>
        tidyr::fill(Quan_Area, .direction = "downup") |>
        dplyr::mutate(QuanMZ = ifelse(Isotope_Label_Type == "Quan", Chromatogram_Precursor_MZ, NA)) |>
        tidyr::fill(QuanMZ, .direction = "downup") |>
        dplyr::mutate(QuanQualRatio = ifelse(Isotope_Label_Type == "Qual", Quan_Area/Area, 1)) |>
        tidyr::replace_na(list(QuanQualRatio = 0)) |>
        dplyr::mutate(QuanQualMZ = paste0(QuanMZ,"/",Chromatogram_Precursor_MZ)) |>

        dplyr::mutate(Quan_Rel_Ab = ifelse(Isotope_Label_Type == "Quan", Rel_Ab, NA)) |>
        tidyr::fill(Quan_Rel_Ab, .direction = "downup") |>
        dplyr::mutate(QuanQual_Rel_Ab_Ratio = ifelse(Isotope_Label_Type == "Qual", Quan_Rel_Ab/Rel_Ab, 1)) |>
        tidyr::replace_na(list(QuanQual_Rel_Ab_Ratio = 0)) |>
        dplyr::ungroup() |>
        dplyr::mutate(MeasVSTheo = QuanQualRatio/QuanQual_Rel_Ab_Ratio) |>
        dplyr::mutate(Is_Outlier = MeasVSTheo > 3 | MeasVSTheo < 0.3) |>
        dplyr::mutate(Is_Outlier = factor(Is_Outlier, levels = c(FALSE, TRUE), labels = c("Within Limit", "Outlier"))) |>
        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, QuanQualMZ, QuanQualRatio, QuanQual_Rel_Ab_Ratio, MeasVSTheo, Is_Outlier) |>
        plotly::plot_ly(x = ~Replicate_Name, y = ~MeasVSTheo,
                        type = 'scatter', mode = 'markers',
                        color = ~Is_Outlier,
                        colors = c('blue', 'red'),
                        text = ~paste("Replicate:", Replicate_Name,
                                      "<br>Homologue Group: ", Molecule,
                                      "<br>Measured against Theoretical Ratio:", round(MeasVSTheo, 1)),
                        marker = list(size = 10)) |>
        layout(title = "Measured/Theoretical ratio >3 or <0.3 are marked in red (ratio of 1 means perfect match",
               xaxis = list(title = "Sample Name"),
               yaxis = list(title = "MeasVSTheo"))

}


##############################################################################

plot_sample_contribution <- function(deconvolution) {

    # How much contribution of each sample to the final deconvoluted homologue group pattern
    plot_data <- deconvolution |>
        tidyr::unnest(deconv_coef) |>
        tidyr::unnest_longer(c(deconv_coef, Batch_Name)) |>
        dplyr::select(Replicate_Name, Batch_Name, deconv_coef)

    # Create the plotly stacked bar plot
    plotly::plot_ly(plot_data,
                    x = ~Replicate_Name,
                    y = ~deconv_coef,
                    type = "bar",
                    color = ~Batch_Name,
                    colors = "Spectral") |>
        plotly::layout(
            title = list(
                text = "Contributions from standards to deconvoluted homologue pattern",
                x = 0.5,  # Center the title
                y = 0.95  # Position slightly down from top
            ),
            barmode = "stack",
            xaxis = list(title = "Replicate Name"),
            yaxis = list(title = "Relative Contribution",
                         tickformat = ".2%"),
            showlegend = TRUE,
            legend = list(title = list(text = "Batch Name"))
        )
}

#############################################################################

plot_homologue_group_pattern_comparison <- function(Sample_distribution, input_selectedSamples){

    # Filter data for selected samples and reshape data
    selected_samples <- Sample_distribution |>
        dplyr::filter(Replicate_Name %in% input_selectedSamples) |>
        dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

    # Get unique homologue groups for consistent coloring
    homologue_groups <- unique(selected_samples$C_homologue)

    # Create a list of plots, one for each Replicate_Name
    plot_list <- selected_samples |>
        split(selected_samples$Replicate_Name) |>
        map(function(df) {
            # Create base plot
            p <- plotly::plot_ly()

            # Add bars for each homologue group
            for(hg in homologue_groups) {
                df_filtered <- df[df$C_homologue == hg,]
                p <- p |>
                    plotly::add_trace(
                        data = df_filtered,
                        x = ~Molecule,
                        y = ~Relative_Area,
                        name = hg,
                        legendgroup = hg,
                        showlegend = (df_filtered$Replicate_Name[1] == input_selectedSamples[1]),
                        type = 'bar',
                        opacity = 1
                    )
            }

            # Add the black line for resolved_distribution
            p <- p |>
                plotly::add_trace(
                    data = df,
                    x = ~Molecule,
                    y = ~resolved_distribution,
                    name = "Deconvoluted Distribution",
                    legendgroup = "DeconvDistr",
                    showlegend = (df$Replicate_Name[1] == input_selectedSamples[1]),
                    type = 'scatter',
                    mode = 'lines+markers',
                    line = list(color = 'black'),
                    marker = list(color = 'black', size = 6),
                    opacity = 0.7
                )

            # Add layout
            p <- p |>
                plotly::layout(
                    xaxis = list(
                        title = "Homologue",
                        tickangle = 45
                    ),
                    yaxis = list(title = "Value"),
                    barmode = 'group',
                    annotations = list(
                        x = 0.5,
                        y = 1.1,
                        text = unique(df$Replicate_Name),
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                    )
                )

            return(p)
        })

    # Combine the plots using subplot
    plotly::subplot(plot_list,
                    nrows = ceiling(length(plot_list)/2),
                    shareX = TRUE,
                    shareY = TRUE) |>
        plotly::layout(
            #title = "Sample Comparison",
            showlegend = TRUE,
            hovermode = 'closest',
            hoverlabel = list(bgcolor = "white"),
            barmode = 'group'
        ) |>
        plotly::config(displayModeBar = TRUE) |>
        htmlwidgets::onRender("
                function(el) {
                    var plotDiv = document.getElementById(el.id);
                    plotDiv.on('plotly_legendclick', function(data) {
                        Plotly.restyle(plotDiv, {
                            visible: data.data[data.curveNumber].visible === 'legendonly' ? true : 'legendonly'
                        }, data.fullData.map((trace, i) => i).filter(i =>
                            data.fullData[i].legendgroup === data.fullData[data.curveNumber].legendgroup
                        ));
                        return false;
                    });
                }
            ")
}



################################################################################
##------------------------------ CPquant.R  ----------------------------------##
# shiny::includeMarkdown("./instructions_CPquant.md")
################################################################################
options(shiny.maxRequestSize = 500 * 1024^2)

ui <- shiny::navbarPage("CPquant",
                        shiny::tabPanel("Quantification Inputs",
                                        shiny::fluidPage(shiny::sidebarLayout(
                                            shiny::sidebarPanel(
                                                width = 3,
                                                shiny::fileInput("fileInput", "Import excel file from Skyline",
                                                                 accept = c('xlsx')),
                                                shiny::textInput("quanUnit", "(Optional) Concentration unit:"),
                                                shiny::radioButtons("quanSum",
                                                                    label = "Choose ions for quantification",
                                                                    choices = c("Quan only", "Sum Quan+Qual"),
                                                                    selected = "Quan only"),
                                                shiny::radioButtons("blankSubtraction",
                                                                    label = "Subtraction with blank?",
                                                                    choices = c("Yes, by avg area of blanks", "No"),
                                                                    selected = "No"),
                                                shiny::radioButtons("correctWithRS",
                                                                    label = "Correct with RS area?",
                                                                    choices = c("Yes", "No"),
                                                                    selected = "No"),
                                                shiny::uiOutput("correctionUI"), # render UI if correctwithRS == "Yes"
                                                shiny::radioButtons("calculateRecovery",
                                                                    label = "Calculate recovery? (req QC samples)",
                                                                    choices = c("Yes", "No"),
                                                                    selected = "No"),
                                                shiny::uiOutput("calculateRecoveryUI"), # render UI if CalculateRecovery == "Yes"
                                                shiny::radioButtons("calculateMDL",
                                                                    label = "Calculate MDL? (req blank samples)",
                                                                    choices = c("Yes", "No"),
                                                                    selected = "No"),
                                                shiny::radioButtons("standardTypes",
                                                                    label = "Types of standards",
                                                                    choices = c("Group Mixtures"),
                                                                    selected = "Group Mixtures"), #work for both single chain and multiple chain mixtures
                                                shiny::tags$div(
                                                    title = "Wait until import file is fully loaded before pressing!",
                                                    shiny::actionButton('go', 'Proceed', width = "100%")
                                                ),
                                                shiny::uiOutput("defineVariables")
                                            ),
                                            shiny::mainPanel(
                                                width = 9,
                                                plotly::plotlyOutput("plot_Skyline_output"),
                                                DT::DTOutput("table_Skyline_output")

                                            )
                                        )
                                        )),
                        shiny::tabPanel(
                            "Input summary",
                            shiny::sidebarPanel(
                                width = 2, # max 12
                                shiny::radioButtons("navSummary", "Choose tab:",
                                                    choices = c("Std Calibration Curves",
                                                                "Removed from Calibration",
                                                                "Quan to Qual ratio",
                                                                "Measured vs Theor Quan/Qual ratio"),
                                                    selected = "Std Calibration Curves"),

                                tags$hr(),
                                tags$p("'Quan to Qual ratio' and 'Measured vs Theor Quan/Qual ratio' only works for Quan only option", style = "font-style: italic;")

                            ),
                            shiny::mainPanel(
                                width = 10,
                                shiny::conditionalPanel(
                                    condition = "input.navSummary == 'Std Calibration Curves'",
                                    tags$h3("Standard calibration curves"),
                                    plotly::plotlyOutput("CalibrationCurves")
                                ),
                                shiny::conditionalPanel(
                                    condition = "input.navSummary == 'Removed from Calibration'",
                                    tags$h3("Calibration series removed from quantification"),
                                    DT::DTOutput("CalibrationRemoved")
                                ),

                                shiny::conditionalPanel(
                                    condition = "input.navSummary == 'Quan to Qual ratio'",
                                    tags$h3("Violin plots of Quant/Qual ions"),
                                    plotly::plotlyOutput("RatioQuantToQual", height = "80vh", width = "100%")
                                ),
                                shiny::conditionalPanel(
                                    condition = "input.navSummary == 'Measured vs Theor Quan/Qual ratio'",
                                    tags$h3("Measured divided by Theoretical Quant/Qual ratios"),
                                    plotly::plotlyOutput("MeasVSTheor", height = "80vh", width = "100%")
                                )
                            )
                        ),

                        shiny::tabPanel(
                            "Quantification summary",
                            shiny::fluidPage(
                                downloadButton("downloadResults", "Export all results to Excel"),
                                shiny::tags$br(), shiny::tags$br(),
                                DT::DTOutput("quantTable"),
                                shiny::tags$br(),
                                plotly::plotlyOutput("sampleContributionPlot")
                            )
                        ),
                        shiny::tabPanel(
                            "Homologue Group Patterns",
                            fluidRow(
                                column(
                                    width = 2,
                                    shiny::radioButtons("plotHomologueGroups", "Choose tab:",
                                                        choices = c("All Samples Overview", "Samples Overlay", "Samples Panels"),
                                                        selected = "All Samples Overview"),
                                    shiny::conditionalPanel(
                                        condition = "input.plotHomologueGroups == 'Samples Overlay'",
                                        shiny::uiOutput("sampleSelectionUIOverlay")
                                    ),
                                    shiny::conditionalPanel(
                                        condition = "input.plotHomologueGroups == 'Samples Panels'",
                                        shiny::uiOutput("sampleSelectionUIComparisons")
                                    ),
                                    shiny::tags$div(
                                        title = "WAIT after pressing..Might take some time before plot shows!",
                                        shiny::actionButton('go2', 'Plot', width = "100%")
                                    )
                                ),
                                column(
                                    width = 10,
                                    # Add a div with margin/padding for spacing
                                    shiny::tags$div(
                                        style = "margin-top: 20px;", # This adds space at the top

                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'All Samples Overview'",
                                            shiny::plotOutput("plotHomologuePatternStatic", height = "80vh", width = "100%")
                                        ),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Overlay'",
                                            plotly::plotlyOutput("plotHomologuePatternOverlay", height = "80vh", width = "100%")
                                        ),
                                        shiny::conditionalPanel(
                                            condition = "input.plotHomologueGroups == 'Samples Panels'",
                                            plotly::plotlyOutput("plotHomologuePatternComparisons", height = "80vh", width = "100%")
                                        )
                                    )
                                )
                            )
                        ),
                        shiny::tabPanel(
                            "QA/QC",
                            shiny::mainPanel(
                                DT::DTOutput("table_recovery"),   # First table output (Skyline recovery data)
                                br(),                     # Optional line break to add space between the tables
                                DT::DTOutput("table_MDL")       # Second table output (LOD table)
                            )
                        ),
                        shiny::tabPanel(
                            "Instructions",
                            shiny::sidebarLayout(
                                shiny::sidebarPanel(shiny::h3("Manual"),
                                                    width = 3),
                                shiny::mainPanel(
                                    shiny::includeMarkdown("./instructions_CPquant.md")
                                )
                            )
                        )

)

################################################################################
server <- function(input, output, session) {

    # define RS for Area correction
    output$correctionUI <- shiny::renderUI({
        if(input$correctWithRS == "Yes") {
            defineCorrectionUI(Skyline_output()) }
    })

    # define RS for recovery
    output$calculateRecoveryUI <- shiny::renderUI({
        if(input$calculateRecovery == "Yes") {
            defineCalcrecoveryUI(Skyline_output()) }
    })

    # Create a reactive object that depends on the quanUnit input
    quantUnit <- reactive({
        input$quanUnit
    })

    Skyline_output <- reactive({
        req(input$fileInput) #requires that the input is available

        # Create a Progress object
        progress <- shiny::Progress$new()
        on.exit(progress$close())

        progress$set(message = "WAIT! Loading data...", value = 0)

        # Read the Skyline output Excel file
        progress$set(value = 0.3, detail = "Reading Excel file")
        df <- readxl::read_excel(input$fileInput$datapath, guess_max = 5000, na = c("", "NA", "#N/A", "N/A")) #guard against different na annotation from Skyline export

        progress$set(value = 0.6, detail = "Processing data")
        # Tidy the input file
        df <- df |>
            dplyr::rename(Replicate_Name = tidyr::any_of(c("Replicate Name", "ReplicateName"))) |> #old versions of Skyline uses space
            dplyr::rename(Sample_Type = tidyr::any_of(c("Sample Type", "SampleType"))) |>
            dplyr::rename(Molecule_List = tidyr::any_of(c("Molecule List", "MoleculeList"))) |>
            dplyr::rename(Mass_Error_PPM = tidyr::any_of(c("Mass Error PPM", "MassErrorPPM"))) |>
            dplyr::rename(Isotope_Label_Type = tidyr::any_of(c("Isotope Label Type", "IsotopeLabelType"))) |>
            dplyr::rename(Chromatogram_Precursor_MZ = tidyr::any_of(c("Chromatogram Precursor M/Z", "ChromatogramPrecursorMz"))) |>
            dplyr::rename(Analyte_Concentration = tidyr::any_of(c("Analyte Concentration", "AnalyteConcentration"))) |>
            dplyr::rename(Batch_Name = tidyr::any_of(c("Batch Name", "BatchName"))) |>
            dplyr::rename(Transition_Note = tidyr::any_of(c("Transition Note", "TransitionNote"))) |>
            dplyr::rename(Sample_Dilution_Factor = tidyr::any_of(c("Sample Dilution Factor", "SampleDilutionFactor"))) |>
            dplyr::mutate(Analyte_Concentration = as.numeric(Analyte_Concentration)) |>
            dplyr::mutate(Area = as.numeric(Area)) |>
            dplyr::mutate(Area = replace_na(Area, 0)) |>
            dplyr::mutate(C_homologue = stringr::str_extract(Molecule, "C\\d+"),
                          Cl_homologue = stringr::str_extract(Molecule, "Cl\\d+"),
                          C_number = as.numeric(stringr::str_extract(C_homologue, "\\d+")),
                          Cl_number = as.numeric(stringr::str_extract(Cl_homologue, "\\d+")),
                          PCA = stringr::str_c(C_homologue, Cl_homologue, sep = "")) |>
            dplyr::mutate(Rel_Ab = as.numeric(map_chr( #Extract relative abundance from second set of curly braces
                Transition_Note, ~ {matches <- str_match_all(.x, "\\{([^}]*)\\}")[[1]]
                if (nrow(matches) >= 2) matches[2, 2] else NA_character_
                })))

        progress$set(value = 0.8, detail = "Applying corrections")

        #EXPERIMENTAL-START
        # Sum all area if quanSum == "Sum Quan+Qual"
        if (input$quanSum == "Sum Quan+Qual") {
            df <- df |>
                dplyr::group_by(Replicate_Name, Molecule_List, Molecule) |>
                dplyr::mutate(Area = sum(Area)) |> #sums Quan and all Qual ions (all ions for the same Molecule will now have same Area)
                dplyr::ungroup()
        }
        #EXPERIMENTAL-END

        # Normalize data based on 'Correct with RS' input
        if (input$correctWithRS == "Yes" && any(df$Molecule_List == "RS")) {
            # Only proceed with RS correction if chooseRS input is available
            if (!is.null(input$chooseRS) && input$chooseRS != "") {
                # Use input$chooseRS directly instead of the reactive
                df <- df |>
                    dplyr::group_by(Replicate_Name) |>
                    dplyr::mutate(Area = Area / first(Area[Molecule == input$chooseRS &
                                                               Molecule_List == "RS" &
                                                               Isotope_Label_Type == "Quan"])) |>
                    dplyr::ungroup()
            }
        }

        # Calculate the average blank value
        if (input$blankSubtraction == "Yes, by avg area of blanks"){
            df_blank <- df |>
                dplyr::filter(Sample_Type == "Blank") |>
                dplyr::group_by(Molecule, Molecule_List, Isotope_Label_Type) |>
                dplyr::summarize(AverageBlank = mean(Area, na.rm = TRUE)) |>
                dplyr::ungroup() |>
                dplyr::filter(!Molecule_List %in% c("IS", "RS", "VS"))

            df <- df |>
                dplyr::full_join(df_blank) |>
                dplyr::mutate(AverageBlank = tidyr::replace_na(AverageBlank, 0)) |>
                dplyr::mutate(Area = dplyr::case_when(Sample_Type == "Unknown" ~ Area - AverageBlank, .default = Area)) |> #only blank subtraction of Unknown Sample Type
                dplyr::mutate(Area = ifelse(Area <0, 0, Area)) #replace negative Area with 0 after blank subtraction
        }


        if (input$standardTypes == "Group Mixtures") {
            df <- df |>
                dplyr::mutate(Quantification_Group = stringr::str_extract(Batch_Name, "^[^_]+")) #extract the initial sequence of characters up to, but not including, the first underscore.
        }

        progress$set(value = 1, detail = "Complete")

        return(df)
    })

    # defineVariablesUI in separate file UI_components.R
    output$defineVariables <- shiny::renderUI({
        defineVariablesUI(Skyline_output())
    })


    # Set reactive values from user input

    removeRsquared <- shiny::eventReactive(input$go, {as.numeric(input$removeRsquared)})
    removeSamples <- shiny::eventReactive(input$go, {as.character(input$removeSamples)})
    Samples_Concentration <- reactiveVal() # Create a reactive value to store deconvolution object into Samples_Concentration() to allow other to access after observeEvent.


    #Render raw table
    output$table_Skyline_output <- DT::renderDT({
        DT::datatable(Skyline_output(),
                      options = list(
                          paging = TRUE,
                          pageLength = 50
                      )
        )
    })


    # Render reactive summary statistics and plots of raw input BEFORE quantification
    # plots.R function
    output$plot_Skyline_output <- plotly::renderPlotly({
        plot_skyline_output(Skyline_output())

    })



    #----------------START: Deconvolution script------------------#

    shiny::observeEvent(input$go, {

        progress <- shiny::Progress$new()
        on.exit(progress$close())

        progress$set(message = "WAIT! Processing data...", value = 0)

        # remove samples if selected by removeSamples input

        if(!is.null(removeSamples()) && length(removeSamples()) > 0){
            Skyline_output_filt <- Skyline_output() |>
                dplyr::filter(!Replicate_Name %in% removeSamples())
        } else{
            Skyline_output_filt <- Skyline_output()
        }



        ##### PREPARE FOR DECONVOLUTION #######
        # Prepare for deconvolution for standards
        progress$set(value = 0.2, detail = "Preparing standards data")

        # Prepare data frame for all standards to be used in deconvolution
        if(input$standardTypes == "Group Mixtures"){

            CPs_standards <- Skyline_output_filt |>
                dplyr::filter(Sample_Type == "Standard",
                              !Molecule_List %in% c("IS", "RS", "VS"), # dont include IS, RS, VS
                              Isotope_Label_Type == "Quan", # use only Quan ions
                              Batch_Name != "NA") |>
                dplyr::mutate(
                    C_range = stringr::str_extract_all(Quantification_Group, "\\d+"), #extract all sequences of digits
                    C_min = as.numeric(purrr::map_chr(C_range, ~.x[1])),
                    C_max = as.numeric(purrr::map_chr(C_range, function(x) {if(length(x) > 1) {x[2]} else {x[1]}}))) |>
                dplyr::select(-C_range) |>
                dplyr::group_by(Batch_Name, Sample_Type, Molecule, Molecule_List, C_number, Cl_number, PCA, Quantification_Group, C_min, C_max) |>
                tidyr::nest() |>
                dplyr::filter(C_number >= C_min & C_number <= C_max) |> # make sure C_number stays within the Quantification_Group chain length
                dplyr::mutate(models = purrr::map(data, ~lm(Area ~ Analyte_Concentration, data = .x))) |> #this includes intercept, if omitting intercept: Area~Analyte_Concentration -1
                dplyr::mutate(coef = purrr::map(models, coef)) |>
                dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"]))|> #get the slope which will be the RF
                dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                dplyr::mutate(cal_rsquared = purrr::map(models, summary)) |> #first create a data frame list with the model
                dplyr::mutate(cal_rsquared = purrr::map(cal_rsquared, purrr::pluck("r.squared"))) |> # then pluck only the r.squared value
                dplyr::select(-coef) |>  # remove coef variable since it has already been plucked
                tidyr::unnest(c(RF, intercept, cal_rsquared)) |>  #removing the list type for these variables
                dplyr::mutate(RF = if_else(RF < 0, 0, RF)) |> # replace negative RF with 0
                dplyr::mutate(cal_rsquared = ifelse(is.nan(cal_rsquared), 0, cal_rsquared)) |>
                dplyr::mutate(RF = if_else(cal_rsquared < removeRsquared(), 0, RF)) |> #replace RF with 0 if rsquared is below removeRsquared()
                dplyr::ungroup() |>
                dplyr::group_by(Batch_Name) |> #grouping by the standards
                dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |> #the sum RF per standard
                dplyr::ungroup()


            # Prepare for deconvolution of samples
            progress$set(value = 0.6, detail = "Preparing sample data")


            CPs_samples <- Skyline_output_filt |>
                dplyr::filter(
                    Sample_Type %in% c("Unknown", "Blank"), #include both unknown and blank
                    !Molecule_List %in% c("IS", "RS", "VS"), # remove IS, RS, VS
                    Isotope_Label_Type == "Quan") |>
                dplyr::group_by(Replicate_Name) |>  # Group by Replicate Name
                dplyr::mutate(Relative_Area = Area / sum(Area, na.rm = TRUE)) |> #Relative area
                dplyr::ungroup() |>
                dplyr::select(-Mass_Error_PPM, -Isotope_Label_Type, -Chromatogram_Precursor_MZ, -Analyte_Concentration, -Batch_Name) |>
                dplyr::mutate(dplyr::across(Relative_Area, ~replace(., is.nan(.), 0)))  # Replace NaN with zero


            CPs_samples_input <- CPs_samples |>
                dplyr::select(Molecule, Replicate_Name, Relative_Area) |>
                tidyr::pivot_wider(names_from = "Replicate_Name", values_from = "Relative_Area")



            # Ensure combined_sample is correctly defined with nested data frames prior to deconvolution
            combined_sample <- CPs_samples  |>
                dplyr::group_by(Replicate_Name, Sample_Type) |>
                tidyr::nest() |>
                dplyr::ungroup()


            ###### Plot calibration curves ######
            # plots.R function
            output$CalibrationCurves <- plotly::renderPlotly({
                plot_calibration_curves(CPs_standards, quantUnit())
            })

            ###### Plot CalibrationRemoved ######
            output$CalibrationRemoved <- DT::renderDT({
                CPs_standards |>
                    dplyr::filter(RF <= 0) |>
                    dplyr::mutate(coef = purrr::map(models, coef)) |>
                    dplyr::mutate(RF = purrr::map_dbl(models, ~ coef(.x)["Analyte_Concentration"]))|> #get the slope which will be the RF
                    dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
                    dplyr::select(Batch_Name, Molecule, Quantification_Group,RF, intercept,  cal_rsquared) |>
                    tidyr::unnest(c(RF, intercept)) |>
                    mutate(across(where(is.numeric), ~ signif(.x, digits = 4))) |>
                    DT::datatable(options = list(pageLength = 40))
            })


            ###### Plot Quan/Qual ratios ######
            # CPquant_plots.R function
            output$RatioQuantToQual <- plotly::renderPlotly({
                plot_quanqualratio(Skyline_output_filt)
            })

            ##### Plot MeasVSTheor ratios ######
            #CPquant_plots.R function
            output$MeasVSTheor <- plotly::renderPlotly({
                plot_meas_vs_theor_ratio(Skyline_output_filt)
            })


            #### DECONVOLUTION ####

            progress$set(value = 0.8, detail = "Performing deconvolution")


            CPs_standards_input <- CPs_standards |>
                dplyr::select(Molecule, Batch_Name, RF) |>
                tidyr::pivot_wider(names_from = Batch_Name, values_from = "RF") |>
                dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))


            CPs_standards_sum_RF <- CPs_standards |>
                dplyr::select(Batch_Name, Sum_RF_group) |>
                dplyr::distinct() |>
                dplyr::ungroup() |>
                tidyr::pivot_wider(names_from = Batch_Name, values_from = "Sum_RF_group") |>
                dplyr::mutate(across(everything(), ~ tidyr::replace_na(., 0)))


            # First populate combined_standard with CPs_standards_input
            combined_standard <- CPs_standards_input  |>
                tibble::column_to_rownames(var = "Molecule") |>
                as.matrix()



            # This performs deconvolution on all mixtures together
            deconvolution <- combined_sample |>
                #perform_deconvolution on only Relative_Area in the nested data frame
                dplyr::mutate(result = purrr::map(data, ~ perform_deconvolution(dplyr::select(.x, Relative_Area), combined_standard, CPs_standards_sum_RF))) |>
                dplyr::mutate(sum_Area = purrr::map_dbl(data, ~sum(.x$Area))) |>
                dplyr::mutate(sum_deconv_RF = as.numeric(purrr::map(result, purrr::pluck("sum_deconv_RF")))) |>
                #EXPERIMENTAL-START
                dplyr::mutate(Sample_Dilution_Factor = purrr::map_dbl(data, ~first(.x$Sample_Dilution_Factor))) |> #since all dilution factor is same for a replicate then take the first
                dplyr::mutate(Concentration = sum_Area/sum_deconv_RF*Sample_Dilution_Factor) |>
                #EXPERIMENTAL-END
                dplyr::mutate(Unit = quantUnit()) |>
                dplyr::mutate(deconv_coef = purrr::map(result, ~as_tibble(list(deconv_coef = .x$deconv_coef, Batch_Name = names(.x$deconv_coef))))) |>
                dplyr::mutate(deconv_rsquared = as.numeric(purrr::map(result, purrr::pluck("deconv_rsquared")))) |>
                dplyr::mutate(deconv_resolved = purrr::map(result, ~tibble::as_tibble(list(deconv_resolved = .x$deconv_resolved, Molecule = rownames(.x$deconv_resolved))))) |>
                dplyr::select(-result)



            #### Calculate the concentration ####

            progress$set(value = 0.9, detail = "Calculating final results")



            # Store Samples_Concentration in the reactive value
            Samples_Concentration(deconvolution)

        }

        progress$set(value = 1, detail = "Complete")


        ### END: Deconvolution script


        # download results from deconvolution to different excel sheets
        output$downloadResults <- shiny::downloadHandler(
            filename = function() {
                paste("CPquant_Results_", Sys.Date(), ".xlsx", sep = "")
            },
            content = function(file) {
                # Create a new workbook
                wb <- openxlsx::createWorkbook()

                # Add worksheets
                openxlsx::addWorksheet(wb, "Quantification")
                openxlsx::addWorksheet(wb, "StandardsContribution")
                openxlsx::addWorksheet(wb, "HomologueDistribution")

                # Write data to worksheets
                openxlsx::writeData(wb, "Quantification",
                                    deconvolution |>
                                        dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |>
                                        dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)))
                openxlsx::writeData(wb, "StandardsContribution",
                                    deconvolution |>
                                        tidyr::unnest(deconv_coef) |>
                                        tidyr::unnest_longer(c(deconv_coef, Batch_Name)) |>
                                        dplyr::select(Replicate_Name, Batch_Name, deconv_coef) |>
                                        dplyr::mutate(deconv_coef = deconv_coef * 100) |>
                                        tidyr::pivot_wider(names_from = Batch_Name, values_from = deconv_coef))
                openxlsx::writeData(wb, "HomologueDistribution",
                                    deconvolution |>
                                        dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                                        dplyr::mutate(data = purrr::map(data, ~ .x |> dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                                        dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                                        tidyr::unnest(data) |>
                                        dplyr::mutate(Deconvoluted_Distribution = as.numeric(resolved_distribution)) |>
                                        dplyr::rename(Relative_Distribution = Relative_Area) |>
                                        dplyr::mutate(Molecule_Concentration = Deconvoluted_Distribution * Concentration) |>
                                        dplyr::select(Replicate_Name, Sample_Type, Molecule_List, Molecule, C_homologue, Cl_homologue, PCA,
                                                      Relative_Distribution, Deconvoluted_Distribution, Molecule_Concentration, Unit)
                )


                # Save workbook
                openxlsx::saveWorkbook(wb, file)
            }
        )


        # Render table
        output$quantTable <- DT::renderDT({
            deconvolution |>
                dplyr::select(Replicate_Name, Sample_Type, Concentration, Unit, deconv_rsquared) |> #select to make compact df for pivot_wider
                dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3)) |>
                DT::datatable(
                    filter = "top", extensions = c("Buttons", "Scroller"),
                    options = list(scrollY = 650,
                                   scrollX = 500,
                                   deferRender = TRUE,
                                   scroller = TRUE,
                                   buttons = list(list(extend = "excel", filename = "Samples_concentration", title = NULL,
                                                       exportOptions = list(
                                                           modifier = list(page = "all")
                                                       )),
                                                  list(extend = "csv", filename = "Samples_concentration", title = NULL,
                                                       exportOptions = list(
                                                           modifier = list(page = "all")
                                                       )),
                                                  list(extend = "colvis", targets = 0, visible = FALSE)),
                                   dom = "lBfrtip",
                                   fixedColumns = TRUE),
                    rownames = FALSE)

        })

        # render sampleContributionPlot
        output$sampleContributionPlot <- renderPlotly({
            plot_sample_contribution(deconvolution)})



        ###########################################################QA/QC#######################################################

        QAQC <- deconvolution |>
            dplyr::select(Replicate_Name, Sample_Type, deconv_rsquared) |>
            dplyr::mutate(deconv_rsquared = round(deconv_rsquared, 3))

        if(input$calculateRecovery == "Yes") {
            # Recovery calculations
            recovery_data <- Skyline_output_filt |>
                dplyr::filter(Isotope_Label_Type == "Quan",
                              Molecule_List %in% c("IS", "RS"))


            RECOVERY <- recovery_data |>  # Calculate recovery
                dplyr::filter(Molecule_List %in% c("RS", "IS") & !(Molecule_List == "RS" & Molecule != input$chooseRS)) |> #remove not chosen RS
                tidyr::pivot_wider(
                    id_cols = c(Replicate_Name, Sample_Type),
                    names_from = Molecule_List,
                    values_from = Area) |>
                dplyr::mutate(across(c(IS, RS), ~replace_na(.x, 0)))

            # Calculate QC ratio
            qc_ratio <- RECOVERY |>
                dplyr::filter(Sample_Type == "Quality Control") |>
                dplyr::mutate(RatioStd = IS / RS) |>
                dplyr::summarize(AverageRatio = mean(RatioStd, na.rm = TRUE))

            # Calculate sample recovery
            RECOVERY <- RECOVERY |>
                dplyr::filter(Sample_Type %in% c("Unknown", "Blank")) |>
                dplyr::mutate(
                    RatioSample = IS / RS,
                    Recovery = RatioSample / as.numeric(qc_ratio$AverageRatio),
                    RecoveryPercentage = round(Recovery * 100, 0)
                ) |>
                dplyr::select(Replicate_Name, Sample_Type, RecoveryPercentage)

            QAQC <- QAQC |>
                dplyr::left_join(RECOVERY, by = c("Replicate_Name", "Sample_Type"))
        }

        # Render recovery table
        output$table_recovery <- DT::renderDT({
            DT::datatable(QAQC,
                          filter = "top",
                          extensions = c("Buttons", "Scroller"),
                          options = list(
                              scrollY = 650,
                              scrollX = 500,
                              deferRender = TRUE,
                              scroller = TRUE,
                              buttons = list(
                                  list(extend = "excel",
                                       filename = "Samples_recovery",
                                       title = NULL,
                                       exportOptions = list(modifier = list(page = "all"))),
                                  list(extend = "csv",
                                       filename = "Samples_recovery",
                                       title = NULL,
                                       exportOptions = list(modifier = list(page = "all"))),
                                  list(extend = "colvis",
                                       targets = 0,
                                       visible = FALSE)
                              ),
                              dom = "lBfrtip",
                              fixedColumns = TRUE
                          ),
                          rownames = FALSE)
        })



        if(input$calculateMDL == "Yes") {
            #MDL calculations (need to take into account if blank subtraction affect or not)
            if (input$blankSubtraction == "No"){

                MDL_data <- deconvolution |>
                    dplyr::filter(Sample_Type == "Blank") |>
                    dplyr::summarize(
                        MDL_sumPCA = mean(Concentration) + 3 * sd(Concentration, na.rm = TRUE),
                        number_of_blanks = dplyr::n_distinct(Replicate_Name)
                    )
            } else if (input$blankSubtraction == "Yes, by avg area of blanks"){
                MDL_data <- deconvolution |>
                    dplyr::filter(Sample_Type == "Blank") |>
                    dplyr::summarize(
                        MDL_sumPCA = 3 * sd(Concentration, na.rm = TRUE),
                        number_of_blanks = dplyr::n_distinct(Replicate_Name)
                    )

            }
        }



        # Render MDL table
        output$table_MDL <- DT::renderDT({
            if (exists("MDL_data")) {
                DT::datatable(MDL_data,
                              options = list(
                                  pageLength = 10,
                                  dom = 't'
                              ),
                              rownames = FALSE)
            } else {
                NULL
            }
        })








    })



    output$sampleSelectionUIOverlay <- renderUI({
        req(Samples_Concentration())

        # Get unique sample names
        sample_names <- unique(Samples_Concentration()$Replicate_Name)

        selectInput("selectedSamples", "Select samples to compare:",
                    choices = sample_names,
                    multiple = TRUE,
                    selected = NULL)
    })

    output$sampleSelectionUIComparisons <- renderUI({
        req(Samples_Concentration())

        # Get unique sample names
        sample_names <- unique(Samples_Concentration()$Replicate_Name)

        selectInput("selectedSamples", "Select samples to compare:",
                    choices = sample_names,
                    multiple = TRUE,
                    selected = NULL)
    })

    shiny::observeEvent(input$go2, {
        req(Samples_Concentration())  # Make sure the data exists

        withProgress(message = 'Generating plot...', value = 0, {

            Sample_distribution <- Samples_Concentration() |>
                dplyr::mutate(data = purrr::map2(data, deconv_resolved, ~dplyr::inner_join(.x, .y, by = "Molecule"))) |>
                dplyr::mutate(data = purrr::map(data, ~ .x |> dplyr::mutate(resolved_distribution = deconv_resolved / sum(deconv_resolved)))) |>
                dplyr::select(-deconv_resolved, -Sample_Dilution_Factor) |>
                tidyr::unnest(data) |>
                dplyr::mutate(resolved_distribution = as.numeric(resolved_distribution)) |>
                dplyr::mutate(deconv_resolved = as.numeric(deconv_resolved)) |>
                dplyr::mutate(Molecule_concentration = resolved_distribution * Concentration) #calculated conc of each molecule


            incProgress(0.5)

            if (input$plotHomologueGroups == "All Samples Overview") {
                output$plotHomologuePatternStatic <- shiny::renderPlot({
                    ggplot2::ggplot(Sample_distribution, aes(x = Molecule, y = Relative_Area, fill = C_homologue)) +
                        ggplot2::geom_col() +
                        ggplot2::facet_wrap(~Replicate_Name) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(axis.text.x = element_blank()) +
                        ggplot2::labs(title = "Relative Distribution (before deconvolution)",
                                      x = "Homologue",
                                      y = "Relative Distribution")
                })
            } else if (input$plotHomologueGroups == "Samples Overlay") {
                output$plotHomologuePatternOverlay <- plotly::renderPlotly({
                    req(input$selectedSamples)
                    req(Sample_distribution)

                    # Filter for selected samples
                    selected_samples <- Sample_distribution |>
                        dplyr::filter(Replicate_Name %in% input$selectedSamples) |>
                        dplyr::mutate(Molecule = factor(Molecule, levels = unique(Molecule[order(C_number, Cl_number)])))

                    req(nrow(selected_samples) > 0)

                    # Create a basic bar plot
                    p <- plotly::plot_ly(data = selected_samples,
                                         x = ~Molecule,
                                         #y = ~resolved_distribution,
                                         y = ~Relative_Area,
                                         color = ~Replicate_Name,
                                         type = "bar",
                                         text = ~paste(
                                             "Sample:", Replicate_Name,
                                             "<br>Homologue:", Molecule,
                                             "<br>Distribution:", round(Relative_Area, 3),
                                             "<br>C-atoms:", C_homologue
                                         ),
                                         hoverinfo = "text"
                    ) |>
                        layout(
                            #title = "Sample Distribution Overlay",
                            xaxis = list(
                                title = "Homologue",
                                tickangle = 45
                            ),
                            yaxis = list(
                                title = "Distribution"
                            ),
                            barmode = 'group',
                            showlegend = TRUE,
                            height = 600,
                            margin = list(b = 100)  # Add more bottom margin for rotated labels
                        )

                    p
                })
            } else if (input$plotHomologueGroups == "Samples Panels") {
                output$plotHomologuePatternComparisons <- plotly::renderPlotly({
                    req(input$selectedSamples)
                    #plots.R function
                    plot_homologue_group_pattern_comparison(Sample_distribution, input$selectedSamples)
                })
            }



            incProgress(1)
        })
    })



    # Close the app when the session ends
    if(!interactive()) {
        session$onSessionEnded(function() {
            stopApp()
            q("no")
        })
    }

}



# Run the application
shinyApp(ui = ui, server = server)
