## ---------------------------- ##
## This is a modified version of the RSA function. The purpose here
## is to provide robust clustered robust standard errors 
## and to use control variables for at least the full model. 
## ---------------------------- ##

RSA.akmod <- function (formula, data = NULL, center = FALSE, scale = FALSE, 
    na.rm = FALSE, out.rm = FALSE, breakline = FALSE, models = "default", 
    cubic = FALSE, verbose = TRUE, add = "", estimator = "MLR", 
    se = "robust", missing = NA, ..., control.variables = c(), cluster.variable = c()) 
{
	require(RSA)
	require(survey)
	require(lavaan.survey)
## ---------------------------- ##
## The original version of the function excludes control variable functionality
## I am going to implement it at least for the full model for now. 
## ---------------------------- ##
#    if (length(control.variables) > 0) 
#        stop("Control.variables feature not implemented yet!")
## ---------------------------- ##

## Editing this to include the fullcluster model ## 

    validmodels <- c("absdiff", "absunc", "diff", "mean", "additive", 
        "IA", "SQD", "SRRR", "SRR", "RR", "SSQD", "SRSQD", "full", 
        "null", "onlyx", "onlyy", "onlyx2", "onlyy2", "weak", 
        "strong", "fullcluster")
    if (length(models) == 1 & models[1] == "all") {
        models <- validmodels
    }
    if (length(models) == 1 & models[1] == "default") {
        models <- c("additive", "IA", "SQD", "SRRR", "SRR", "RR", 
            "SSQD", "SRSQD", "full", "null", "onlyx2", "onlyy2", 
            "onlyx", "onlyy")
    }
    if (any(!models %in% validmodels)) 
        stop("Unknown model name provided in parameter 'models'.")
    s.NULL <- s.full <- s.full.cluster <- s.IA <- s.diff <- s.mean <- s.absdiff <- s.additive <- s.SQD <- s.SSQD <- s.SRSQD <- s.absunc <- s.cubic <- s.RR <- s.SRR <- s.SRRR <- s.onlyx <- s.onlyy <- s.onlyx2 <- s.onlyy2 <- s.weak <- s.strong <- NULL
    SRSQD.rot <- ""
    SRRR.rot <- ""
    add <- paste0("\n# User defined syntax:\n", add)
    
## ---------------------------- ##
## This section of RSA creates scaled variables, creates the polynomial terms, checks the range of variables
## and checks for missing values 
## ---------------------------- ##    
    DV <- all.vars(formula)[1]
    IV1 <- all.vars(formula)[2]
    IV2 <- all.vars(formula)[3]
    df <- data[, c(DV, IV1, IV2, control.variables, cluster.variable)]
    df[, IV1] <- scale(df[, IV1], center = center, scale = scale)
    df[, IV2] <- scale(df[, IV2], center = center, scale = scale)
    df <- add.variables(formula, data.frame(data.matrix(df)))
    if (0 < min(df[, IV1], na.rm = TRUE) | 0 > max(df[, IV1], 
        na.rm = TRUE)) 
        warning(paste("The numerical zero point is outside of the range of variable", 
            IV1, ". Please consider re-centering the variable."))
    if (0 < min(df[, IV2], na.rm = TRUE) | 0 > max(df[, IV2], 
        na.rm = TRUE)) 
        warning(paste("The numerical zero point is outside of the range of variable", 
            IV2, ". Please consider re-centering the variable."))
    if ((max(df[, IV1], na.rm = TRUE) - min(df[, IV1], na.rm = TRUE))/(max(df[, 
        IV2], na.rm = TRUE) - min(df[, IV2], na.rm = TRUE)) > 
        2) 
        warning("Predictor variables have a very different range (by factor 2 or larger)- please check scaling of variables.")
    if (is.na(missing)) {
        if (any(is.na(df))) {
            missing <- "fiml"
            warning("There are missing values in your data set. Model is computed with option `missing = 'fiml'`. This is only valid if the data are missing completely at random (MCAR) or missing at random (MAR)! If you want to exclude NA, use `missing = 'listwise'`", 
                call. = FALSE)
        }
        else {
            missing <- "listwise"
        }
    }
## ---------------------------- ##
## This section of RSA creates the string names 
## of the newly created variables (above) for higher order terms and 
## interaction terms. This also creates the addition for control variables. 
## ---------------------------- ##    
    IV12 <- paste0(IV1, "2")
    IV22 <- paste0(IV2, "2")
    IV13 <- paste0(IV1, "3")
    IV23 <- paste0(IV2, "3")
    IV_IA <- paste0(IV1, "_", IV2)
    IV_IA2 <- paste0(IV1, "_", IV2, "2")
    IV_IA3 <- paste0(IV1, "2", "_", IV2)
    W_IV1 <- paste0("W_", IV1)
    W_IV2 <- paste0("W_", IV2)
    CV <- ifelse(length(control.variables > 0), paste0(" + ", 
        paste(control.variables, collapse = " + ")), "")
    addcubic <- ""
    if (cubic == TRUE) 
        addcubic <- paste0(" + ", paste(IV13, IV23, IV_IA2, IV_IA3, 
            sep = " + "))
    f <- paste0(paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA, 
        IV22, sep = " + ")), addcubic, CV)
        
## ---------------------------- ##
# This uses regression to get model statistics and examine for outliers
## ---------------------------- ##        

## ---------------------------- ##        
# AK NOTE: Need to modify this to provide the summary statistics for the model
# that has the control variables included. Maybe include an additional model
# that way we can have a change in the F r2 test from control to 
# inclusion of the polynomial terms.
## ---------------------------- ##        

    lm.full <- lm(f, df, na.action = na.exclude)
    if (is.null(out.rm) || (typeof(out.rm) == "logical" && out.rm == 
        TRUE)) {
#        out.rm <- "bj1980"
		out.rm = "none"
    }
    if ((typeof(out.rm) == "logical" && out.rm == FALSE)) {
        out.rm <- "none"
    }
    out.rm <- match.arg(out.rm, c("bj1980", "robust", "none"))
    df$out <- FALSE
    if (out.rm == "bj1980") {
        inf <- influence.measures(lm.full)
        df$out <- apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 
            1, sum) == 3
        n.out <- sum(na.omit(df$out) == TRUE)
        if (verbose == TRUE & n.out > 0) {
            warning(paste("Removed", n.out, "multivariate outlier(s) according to Bollen & Jackman (1980) criteria. Outliers are in row(s):", 
                paste(which(df$out == TRUE), collapse = ", ")))
        }
    }
    if (out.rm == "robust") {
        stop("Robust outlier detection not implemented yet.")
    }
    df$out[is.na(df$out)] <- FALSE
## ---------------------------- ##
# This section of RSA builds the polynomial equations and runs the 
# path analysis. 
## ---------------------------- ##    
    
    withCallingHandlers({
        poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", 
            IV12, " + b4*", IV_IA, " + b5*", IV22, CV)
        if ("null" %in% models) {
            s.NULL <- sem(paste0(DV, "~ 1 + 0*", IV1, " + 0*", 
                IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", 
                IV22, CV), data = df[df$out == FALSE, ], fixed.x = TRUE, 
                meanstructure = TRUE, se = se, estimator = estimator, 
                missing = missing, ...)
        }
        if ("additive" %in% models) {
            if (verbose == TRUE) 
                print("Computing additive model (additive) ...")
            m.additive <- paste(poly, "b3==0", "b4==0", "b5==0", 
                "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", add, sep = "\n")
            s.additive <- sem(m.additive, data = df[df$out == 
                FALSE, ], fixed.x = TRUE, meanstructure = TRUE, 
                se = se, estimator = estimator, missing = missing, 
                ...)
        }
        if ("onlyx2" %in% models) {
            if (verbose == TRUE) 
                print("Computing x + x^2 model (onlyx2) ...")
            m.onlyx2 <- paste(poly, "b2==0", "b4==0", "b5==0", 
                "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", add, sep = "\n")
            s.onlyx2 <- sem(m.onlyx2, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("onlyy2" %in% models) {
            if (verbose == TRUE) 
                print("Computing y + y^2 model (onlyy2) ...")
            m.onlyy2 <- paste(poly, "b1==0", "b3==0", "b4==0", 
                "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", add, sep = "\n")
            s.onlyy2 <- sem(m.onlyy2, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("onlyx" %in% models) {
            if (verbose == TRUE) 
                print("Computing x model (onlyx) ...")
            m.onlyx <- paste(poly, "b2==0", "b3==0", "b4==0", 
                "b5==0", "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", add, sep = "\n")
            s.onlyx <- sem(m.onlyx, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("onlyy" %in% models) {
            if (verbose == TRUE) 
                print("Computing y model (onlyy) ...")
            m.onlyy <- paste(poly, "b1==0", "b3==0", "b4==0", 
                "b5==0", "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", add, sep = "\n")
            s.onlyy <- sem(m.onlyy, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("diff" %in% models) {
            if (verbose == TRUE) 
                print("Computing difference model (diff) ...")
            m.diff <- paste(poly, "b3==0", "b4==0", "b5==0", 
                "b1 == -b2", "a1 := b1+b2", "a2 := 0", "a3 := b1-b2", 
                "a4 := 0", add, sep = "\n")
            s.diff <- sem(m.diff, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("mean" %in% models) {
            if (verbose == TRUE) 
                print("Computing mean model (mean) ...")
            m.mean <- paste(poly, "b3==0", "b4==0", "b5==0", 
                "b1 == b2", "a1 := b1+b2", "a2 := 0", "a3 := b1-b2", 
                "a4 := 0", add, sep = "\n")
            s.mean <- sem(m.mean, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("IA" %in% models) {
            if (verbose == TRUE) 
                print("Computing interaction model (IA)...")
            m.IA <- paste(poly, "b3==0", "b5==0", "a1 := b1+b2", 
                "a2 := b3+b4+b5", "a3 := b1-b2", "a4 := b3-b4+b5", 
                "a5 := b3-b5", "X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)", 
                "Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "p10 := Y0 - p11*X0", "p20 := Y0 - p21*X0", "PA1.curv := b3 + b4*p11 + b5*(p11^2)", 
                "PA2.curv := b3 + b4*p21 + b5*(p21^2)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.IA <- sem(m.IA, data = df[df$out == FALSE, ], fixed.x = TRUE, 
                meanstructure = TRUE, se = se, estimator = estimator, 
                missing = missing, ...)
        }
        if ("SQD" %in% models) {
            if (verbose == TRUE) 
                print("Computing squared difference model (SQD) ...")
            m.SQD <- paste(poly, "b1==0", "b2==0", "b3==b5", 
                "b3+b4+b5==0", "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SQD <- sem(m.SQD, data = df[df$out == FALSE, ], 
                fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("SSQD" %in% models) {
            if (verbose == TRUE) 
                print("Computing shifted squared difference model (SSQD) ...")
            m.SSQD <- paste(poly, "b1==-b2", "b3==b5", "b3+b4+b5==0", 
                "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "C := b1 / (2*b3)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SSQD <- sem(m.SSQD, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if (any(models %in% c("RR"))) {
            if (verbose == TRUE) 
                print("Computing rising ridge model (RR) ...")
            m.RR <- paste(poly, "b1==b2", "b3==b5", "b3+b4+b5==0", 
                "a1 := b1+b2", "a2 := b3+b4+b5", "a3 := b1-b2", 
                "a4 := b3-b4+b5", "a5 := b3-b5", "meaneffect := b1+b2", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.RR <- sem(m.RR, data = df[df$out == FALSE, ], fixed.x = TRUE, 
                meanstructure = TRUE, se = se, estimator = estimator, 
                missing = missing, ...)
        }
        if (any(models %in% c("SRR"))) {
            if (verbose == TRUE) 
                print("Computing shifted rising ridge model (SRR) ...")
            m.SRR <- paste(poly, "b3==b5", "b3+b4+b5==0", "a1 := b1+b2", 
                "a2 := b3+b4+b5", "a3 := b1-b2", "a4 := b3-b4+b5", 
                "a5 := b3-b5", "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "meaneffect := a1", "C := (b1-b2) / (4*b3)", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SRR <- sem(m.SRR, data = df[df$out == FALSE, ], 
                fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if (any(models %in% c("SRRR"))) {
            if (verbose == TRUE) 
                print("Computing rotated and shifted rising ridge model (SRRR), up ...")
            m.SRRR.up <- paste(paste(poly, " + start(0.01)*", 
                IV12, " + start(0.01)*", IV22), "b3 > 0.000001", 
                "b5 > 0.000001", "b4^2 == 4*b3*b5", "a1 := b1+b2", 
                "a2 := b3+b4+b5", "a3 := b1-b2", "a4 := b3-b4+b5", 
                "a5 := b3-b5", "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "meaneffect := (b2*b4 - 2*b1*b5) / b4", "C := (-2*b1*b5 - b2*b4) / (4*b4*b5)", 
                "S := (-b4) / (2*b5)", "a4.rescaled := b3/S^2 - b4/S + b5", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SRRR.up <- sem(m.SRRR.up, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
            if (verbose == TRUE) 
                print("Computing rotated and shifted rising ridge model (SRRR), down ...")
            m.SRRR.down <- paste(paste(poly, " + start(-0.01)*", 
                IV12, " + start(-0.01)*", IV22), "b3 < -0.000001", 
                "b5 < -0.000001", "b4^2 == 4*b3*b5", "a1 := b1+b2", 
                "a2 := b3+b4+b5", "a3 := b1-b2", "a4 := b3-b4+b5", 
                "a5 := b3-b5", "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "meaneffect := (b2*b4 - 2*b1*b5) / b4", "C := (-2*b1*b5 - b2*b4) / (4*b4*b5)", 
                "S := (-b4) / (2*b5)", "a4.rescaled := b3/S^2 - b4/S + b5", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SRRR.down <- sem(m.SRRR.down, data = df[df$out == 
                FALSE, ], fixed.x = TRUE, meanstructure = TRUE, 
                se = se, estimator = estimator, missing = missing, 
                ...)
            if (inspect(s.SRRR.up, "converged") == FALSE & inspect(s.SRRR.down, 
                "converged") == TRUE) {
                SRRR.rot <- "down"
            }
            else if (inspect(s.SRRR.up, "converged") == TRUE & 
                inspect(s.SRRR.down, "converged") == FALSE) {
                SRRR.rot <- "up"
            }
            else if (inspect(s.SRRR.up, "converged") == TRUE & 
                inspect(s.SRRR.down, "converged") == TRUE) {
                SRRR.rot <- ifelse(fitMeasures(s.SRRR.up, "chisq") > 
                  fitMeasures(s.SRRR.down, "chisq"), "down", 
                  "up")
            }
            else {
                if (verbose == TRUE) 
                  print("Warning: SRRR model has not converged (neither up nor down curvature)")
            }
            if (SRRR.rot == "up") {
                s.SRRR <- s.SRRR.up
            }
            else if (SRRR.rot == "down") {
                s.SRRR <- s.SRRR.down
            }
            if (verbose == TRUE) 
                print(paste0("Direction of SRRR curvature: ", 
                  SRRR.rot))
        }
        if (any(models %in% c("SRSQD"))) {
            if (verbose == TRUE) 
                print("Computing rotated squared difference model (SRSQD), up ...")
            m.SRSQD.up <- paste(paste(poly, " + start(0.001)*", 
                IV22), "b1 == (b2*b4)/(2*b5)", "b3 > 0.000001", 
                "b5 > 0.000001", "b4^2 == 4*b3*b5", "C := -.5*(b2/b5)", 
                "S := (-b4) / (2*b5)", "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "a4.rescaled := b3/S^2 - b4/S + b5", "X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)", 
                "Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "p10 := Y0 - p11*X0", "p20 := Y0 - p21*X0", "PA1.curv := b3 + b4*p11 + b5*(p11^2)", 
                "PA2.curv := b3 + b4*p21 + b5*(p21^2)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SRSQD.up <- sem(m.SRSQD.up, data = df[df$out == 
                FALSE, ], fixed.x = TRUE, meanstructure = TRUE, 
                se = se, estimator = estimator, missing = missing, 
                ...)
            if (verbose == TRUE) 
                print("Computing rotated squared difference model (SRSQD), down ...")
            m.SRSQD.down <- paste(paste(poly, " + start(-0.001)*", 
                IV22), "b1 == (b2*b4)/(2*b5)", "b3 < -0.000001", 
                "b5 < -0.000001", "b4^2 == 4*b3*b5", "C := -.5*(b2/b5)", 
                "S := (-b4) / (2*b5)", "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "a4.rescaled := b3/S^2 - b4/S + b5", "X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)", 
                "Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p10 := Y0 - p11*X0", "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "p20 := Y0 - p21*X0", "PA1.curv := b3 + b4*p11 + b5*(p11^2)", 
                "PA2.curv := b3 + b4*p21 + b5*(p21^2)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                add, sep = "\n")
            s.SRSQD.down <- sem(m.SRSQD.down, data = df[df$out == 
                FALSE, ], fixed.x = TRUE, meanstructure = TRUE, 
                se = se, estimator = estimator, missing = missing, 
                ...)
            if (inspect(s.SRSQD.up, "converged") == FALSE & inspect(s.SRSQD.down, 
                "converged") == TRUE) {
                SRSQD.rot <- "down"
            }
            else if (inspect(s.SRSQD.up, "converged") == TRUE & 
                inspect(s.SRSQD.down, "converged") == FALSE) {
                SRSQD.rot <- "up"
            }
            else if (inspect(s.SRSQD.up, "converged") == TRUE & 
                inspect(s.SRSQD.down, "converged") == TRUE) {
                SRSQD.rot <- ifelse(fitMeasures(s.SRSQD.up, "chisq") > 
                  fitMeasures(s.SRSQD.down, "chisq"), "down", 
                  "up")
            }
            else {
                if (verbose == TRUE) 
                  warning("Warning: SRSQD model has not converged (neither up nor down curvature)")
            }
            if (SRSQD.rot == "up") {
                s.SRSQD <- s.SRSQD.up
            }
            else if (SRSQD.rot == "down") {
                s.SRSQD <- s.SRSQD.down
            }
            if (verbose == TRUE) 
                print(paste0("Direction of SRSQD curvature: ", 
                  SRSQD.rot))
        }
## ---------------------------- ##
## Here is the polynomial model that I'm going to alter. It is going to use
## clustered robust standard errors (if the user specified a clustering variable
## ---------------------------- ##        

        if ("full" %in% models) {
            if (verbose == TRUE) 
                print("Computing polynomial model (full) ...")
        	m.full <- paste(poly, "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)", 
                "Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p10 := Y0 - p11*X0", "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "p20 := Y0 - p21*X0", "PA1.curv := b3 + b4*p11 + b5*(p11^2)", 
                "PA2.curv := b3 + b4*p21 + b5*(p21^2)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "weakcondition    := b3*b5", "strongcondition1 := (b2*b4)/(2*b5) - b1", 
                "strongcondition2 := 2*sqrt(b3*b5)  - b4", add, 
                sep = "\n")
                
            # This model is not going to deal with missing values in the way that is specified above. It will just use the default for SEM, which is dependent upon the type of estimator that is used. 
            
            # Have to actually create a full string of this so that the full call is included in the s.full object. If I just use m.full in a regular sem call, then I'll get m.full in the output object
            
             call.full = paste("sem(model='",m.full,"', data=df[df$out == FALSE, ], fixed.x=TRUE, meanstructure=TRUE, se='",se,"', estimator='",estimator,"', ...)", sep="")
             
			 str_eval <- function(x) {return(eval(parse(text=x)))}
			 s.full <- str_eval(call.full)
             			
			## ------------------ ##
			## This is the only change, but it creates an additional model
			## to report
			## ------------------ ##			
			                
            if("fullcluster" %in% models) {				            	
            	d2 = svydesign(ids=~get(cluster.variable), data=df)           	
				s.full.cluster = lavaan.survey(s.full, survey.design=d2, estimator=estimator) 
           	
            }                
        }
        
        if ("weak" %in% models) {
            if (verbose == TRUE) 
                print("Computing weak fit pattern ...")
            m.weak <- paste(poly, "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)", 
                "Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p10 := Y0 - p11*X0", "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "p20 := Y0 - p21*X0", "PA1.curv := b3 + b4*p11 + b5*(p11^2)", 
                "PA2.curv := b3 + b4*p21 + b5*(p21^2)", "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "b3*b5 > 0", add, sep = "\n")
            s.weak <- sem(m.weak, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("strong" %in% models) {
            if (verbose == TRUE) 
                print("Computing strong fit pattern ...")
            m.strong <- paste(poly, "a1 := b1+b2", "a2 := b3+b4+b5", 
                "a3 := b1-b2", "a4 := b3-b4+b5", "a5 := b3-b5", 
                "p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4", 
                "p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
                "PA1.curv := b3 + b4*p11 + b5*(p11^2)", "PA2.curv := b3 + b4*p21 + b5*(p21^2)", 
                "l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
                "b3*b5 > 0.000001", "(b2*b4) == 2*b1*b5", "4*b3*b5  == b4^2", 
                add, sep = "\n")
            s.strong <- sem(m.strong, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if (cubic == TRUE) {
            if (verbose == TRUE) 
                print("Computing full cubic model (cubic) ...")
            m.cubic <- paste(paste0(poly, " + b9*", IV13, " + b10*", 
                IV_IA2, " + b11*", IV_IA3, " + b12*", IV23), 
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b9 + b10 + b11 + b12", 
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b9 + b10 - b11 - b12", 
                add, sep = "\n")
            s.cubic <- sem(m.cubic, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("absdiff" %in% models) {
            if (verbose == TRUE) 
                print("Computing constrained absolute difference model (absdiff) ...")
            m.absdiff <- paste(paste0(DV, " ~ b1*", IV1, " + b2*", 
                IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2), 
                "b1 == 0", "b2 == 0", "b6 == 0", "b7 == -b8", 
                add, sep = "\n")
            s.absdiff <- sem(m.absdiff, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
        if ("absunc" %in% models) {
            if (verbose == TRUE) 
                print("Computing unconstrained absolute difference model (absunc) ...")
            m.absunc <- paste(paste0(DV, " ~ b1*", IV1, " + b2*", 
                IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2), 
                ifelse(breakline == FALSE, "b6==0", ""), add, 
                sep = "\n")
            s.absunc <- sem(m.absunc, data = df[df$out == FALSE, 
                ], fixed.x = TRUE, meanstructure = TRUE, se = se, 
                estimator = estimator, missing = missing, ...)
        }
    }, warning = function(w) {
        W <- as.character(w$call)
        if ((W[1] == "sqrt" & W[2] == "diag(def.cov)" & grepl("NaNs", 
            w$message)) | (W[1] == "sqrt") | (W[1] == "nlminb" & 
            W[2] == "x.par") | (W[2] %in% c("m.SRRR.up", "m.SRRR.down", 
            "m.SRSQD.up", "m.SRSQD.down") & grepl("model has NOT converged", 
            w$message))) {
            invokeRestart("muffleWarning")
        }
    })
    chisq1 <- plyr::ldply(list(full = s.full, SRRR = s.SRRR, 
        SRR = s.SRR, RR = s.RR, SQD = s.SQD), function(x) {
        chi <- -1
        if (!is.null(x)) {
            if (inspect(x, "converged") == TRUE) 
                chi <- fitMeasures(x, "chisq")
        }
        return(chi)
    })

    chisq1 <- chisq1[chisq1[, 2] >= 0, ]
    if (nrow(chisq1) > 1) {
        chisq1$lag <- c(diff(chisq1[, 2], lag = 1), NA)
        if (any(chisq1$lag < 0, na.rm = TRUE)) {
            warning(paste0("There are convergence problems with model ", 
                chisq1[which(chisq1$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
        }
    }
    chisq2 <- plyr::ldply(list(full = s.full, SRRR = s.SRRR, 
        SRSQD = s.SRSQD, SSQD = s.SSQD, SQD = s.SQD), function(x) {
        chi <- -1
        if (!is.null(x)) {
            if (inspect(x, "converged") == TRUE) 
                chi <- fitMeasures(x, "chisq")
        }
        return(chi)
    })
    chisq2 <- chisq2[chisq2[, 2] >= 0, ]
    if (nrow(chisq1) > 1) {
        chisq2$lag <- c(diff(chisq2[, 2], lag = 1), NA)
        if (any(chisq2$lag < 0, na.rm = TRUE)) {
            warning(paste0("There are convergence problems with model ", 
                chisq2[which(chisq2$lag < 0), ".id"], ". Its chi-square value is higher than that of a nested model, which is theoretically not possible. Please inspect the results with care, using the compare()-function"))
        }
    }
    modellist <- list(null = s.NULL, full = s.full, fullcluster = s.full.cluster, IA = s.IA, 
        diff = s.diff, mean = s.mean, absdiff = s.absdiff, additive = s.additive, 
        SQD = s.SQD, SRRR = s.SRRR, SRR = s.SRR, RR = s.RR, SSQD = s.SSQD, 
        SRSQD = s.SRSQD, absunc = s.absunc, cubic = s.cubic, 
        onlyx = s.onlyx, onlyy = s.onlyy, onlyx2 = s.onlyx2, 
        onlyy2 = s.onlyy2, weak = s.weak, strong = s.strong)
    res <- list(models = modellist, SRSQD.rot = SRSQD.rot, SRRR.rot = SRRR.rot, 
        LM = summary(lm.full), formula = formula, data = df, 
        out.rm = out.rm, outliers = which(df$out == TRUE), DV = DV, 
        IV1 = IV1, IV2 = IV2, IV12 = IV12, IV22 = IV22, IV_IA = IV_IA, 
        W_IV1 = W_IV1, W_IV2 = W_IV2, IV13 = IV13, IV23 = IV23, 
        IV_IA2 = IV_IA2, IV_IA3 = IV_IA3, r.squared = summary(lm.full)$r.squared)
    attr(res, "class") <- "RSA"
    return(res)
}
environment(RSA.akmod) <- asNamespace('RSA')