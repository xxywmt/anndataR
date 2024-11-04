#' Convert a Seurat object to an AnnData object
#'
#' `to_Seurat()` converts an AnnData object to a Seurat object.
#'
#' @param obj An AnnData object
#'
#' @importFrom Matrix t
#'
#' @noRd
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' to_Seurat(ad)
# TODO: Add parameters to choose which how X and layers are translated into counts, data and scaled.data
to_Seurat <- function(obj) { # nolint
  requireNamespace("SeuratObject")

  stopifnot(inherits(obj, "AbstractAnnData"))

  # translate var_names
  # trackstatus: class=Seurat, feature=get_var_names, status=done
  var_names_ <- .toseurat_check_obsvar_names(obj$var_names, "var_names")

  # translate obs_names
  # trackstatus: class=Seurat, feature=get_obs_names, status=done
  obs_names_ <- .toseurat_check_obsvar_names(obj$obs_names, "obs_names")

  # translate var
  # trackstatus: class=Seurat, feature=get_var, status=done
  var_ <- obj$var
  rownames(var_) <- var_names_

  # translate obs
  # trackstatus: class=Seurat, feature=get_obs, status=done
  obs_ <-
    if (ncol(obj$obs) > 0) {
      ob <- obj$obs
      rownames(ob) <- obs_names_
      ob
    } else {
      NULL
    }

  # translate X
  # trackstatus: class=Seurat, feature=get_X, status=wip
  # TODO: should x_ be passed to counts or to data?
  # TODO: creating a seurat object when th AnnData doesn't contain X or layers
  # probably doesn't make any sense
  x_ <-
    if (!is.null(obj$X)) {
      Matrix::t(obj$X)
    } else {
      mat <- Matrix::sparseMatrix(
        i = integer(0),
        p = c(0L),
        x = integer(0),
        dims = c(obj$n_vars(), obj$n_obs())
      )
      attr(mat, "is_X_null") <- TRUE # nolint
      mat
    }
  dimnames(x_) <- list(var_names_, obs_names_)
  x_assay <- SeuratObject::CreateAssayObject(counts = x_)

  # create seurat object
  if (ncol(var_) > 0) {
    # don't add var metadata if the data frame does not contain any columns
    x_assay <- SeuratObject::AddMetaData(x_assay, metadata = var_)
  }
  seurat_obj <- SeuratObject::CreateSeuratObject(x_assay, meta.data = obs_)

  # add layers
  # trackstatus: class=Seurat, feature=get_layers, status=wip
  # TODO: should values be passed to counts or to data?
  for (key in obj$layers_keys()) {
    layer_ <- t(obj$layers[[key]])
    dimnames(layer_) <- list(var_names_, obs_names_)
    seurat_obj[[key]] <- SeuratObject::CreateAssayObject(counts = layer_)
  }

  seurat_obj
}

.toseurat_check_obsvar_names <- function(names, label) {
  if (any(grepl("_", names))) {
    # mimic seurat behaviour
    warning(wrap_message(
      "'", label, "' ",
      "cannot have underscores ('_') when converting to Seurat, ",
      "replacing with dashes ('-')"
    ))
    names <- gsub("_", "-", names)
  }

  names
}

#' Convert a Seurat object to an AnnData object
#'
#' `from_Seurat()` converts a Seurat object to an AnnData object.
#' Only one assay can be converted at a time.
#'
#' For more information on the functionality of an AnnData object, see [anndataR-package].
#'
#' @param seurat_obj An object inheriting from Seurat.
#' @param output_class Name of the AnnData class. Must be one of `"HDF5AnnData"` or `"InMemoryAnnData"`.
#' @param assay Assay to be converted. If NULL, `DefaultAssay()` is used.
#' @param X Which of 'counts', 'data', or 'scale.data' will be used for X. By default, 'counts' will be used (if it is
#'   not empty), followed by 'data', then 'scale.data'. The remaining non-empty slots will be stored in different
#'   layers.
#' @param ... Additional arguments passed to the generator function.
#' @param layers
#'
#' @export
#'
#' @seealso [anndataR-package]
# TODO: Add examples
# nolint start: object_name_linter
from_Seurat <- function(
    # nolint end: object_name_linter
    seurat_obj,
    output_class = c("InMemoryAnnData", "HDF5AnnData"),
    assay = NULL,
    X = "counts",
    layers = c("counts", "data", "scale.data"),
    ...) {
  output_class <- match.arg(output_class)

  stopifnot(inherits(seurat_obj, "Seurat"))

  if (!is.null(X)) {
    if (!X %in% c("counts", "data", "scale.data")) {
      stop("X must be NULL or one of: 'counts', 'data', 'scale.data'")
    }
  }

  # If a specific assay is selected, use it
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop("'assay' must be NULL or one of: ", paste0("'", names(seurat_obj@assays), "'", collapse = ", "))
    }
    assay_name <- assay
  } else {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)

    # If Seurat object contains multiple assays, notify user the Default one is used
    if (length(names(seurat_obj@assays)) > 1) {
      message(
        "There are ", length(names(seurat_obj@assays)), " assays in the Seurat object; using the default assay ('",
        assay_name, "'). You can use the `assay` parameter to select a specific assay."
      )
    }
  }

  # get obs
  # trackstatus: class=Seurat, feature=set_obs_names, status=done
  # trackstatus: class=Seurat, feature=set_obs, status=done
  obs <- seurat_obj@meta.data
  rownames(obs) <- colnames(seurat_obj) # TODO: this is probably not needed

  # construct var
  # trackstatus: class=Seurat, feature=set_var_names, status=done
  # trackstatus: class=Seurat, feature=set_var, status=done
  var <- seurat_obj@assays[[assay_name]]@meta.features
  rownames(var) <- rownames(seurat_obj@assays[[assay_name]]) # TODO: this is probably not needed

  # use generator to create new AnnData object
  generator <- get_anndata_constructor(output_class)
  ad <- generator$new(
    obs = obs,
    var = var,
    ...
  )

  if (!is.null(X)) {
    # Check if the slot is not empty
    if (all(dim(SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)) == 0)) {
      stop("The '", X, "' slot is empty.")
    }

    assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = X, assay = assay_name)

    # Remove names
    dimnames(assay_data) <- list(NULL, NULL)
    ad$X <- Matrix::t(assay_data)
  } else {
    # Cannot compare other values with NULL
    X <- "none"
  }

  # Add the remaining non-empty slots as layers
  slots <- layers
  slots <- slots[slots != X]

  for (slot in slots) {
    if (!all(dim(SeuratObject::GetAssayData(seurat_obj, slot = slot)) == 0)) {
      assay_data <- SeuratObject::GetAssayData(seurat_obj, slot = slot, assay = assay_name)
      dimnames(assay_data) <- list(NULL, NULL)
      ad$layers[[slot]] <- Matrix::t(assay_data)
    }
  }

  return(ad)
}

#' Convert a seurat object to an anndata object using reticulate
#'
#' This function takes a Seurat object and converts it to an anndata object using the reticulate package.
#'
#' @param srt A Seurat object.
#' @param assay_X Assay to convert as the main data matrix (X) in the anndata object.
#' @param slot_X Slot name for assay_X in the Seurat object.
#' @param assay_layers Assays to convert as layers in the anndata object.
#' @param slot_layers Slot names for the assay_layers in the Seurat object.
#' @param convert_tools Logical indicating whether to convert the tool-specific data.
#' @param convert_misc Logical indicating whether to convert the miscellaneous data.
#' @param features Optional vector of features to include in the anndata object. Defaults to all features in assay_X.
#' @param verbose Logical indicating whether to print verbose messages during the conversion process.
#'
#' @return An anndata object.
#'
#' @return A \code{anndata} object.
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata
#'
#' ### Or save as an h5ad file or a loom file
#' # adata$write_h5ad("pancreas_sub.h5ad")
#' # adata$write_loom("pancreas_sub.loom", write_obsm_varm = TRUE)
#'
#' @importFrom reticulate import np_array
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix t
#' @export
srt_to_adata <- function(srt, features = NULL,
                         assay_X = "RNA", slot_X = "counts",
                         assay_layers = c("spliced", "unspliced"), slot_layers = "counts",
                         convert_tools = FALSE, convert_misc = FALSE, verbose = TRUE) {
  #check_Python(c("scanpy", "numpy"))

  if (!inherits(srt, "Seurat")) {
    stop("'srt' is not a Seurat object.")
  }

  if (is.null(features)) {
    features <- rownames(srt[[assay_X]])
  }
  if (length(slot_layers) == 1) {
    slot_layers <- rep(slot_layers, length(assay_layers))
    names(slot_layers) <- assay_layers
  } else if (length(slot_layers) != length(assay_layers)) {
    stop("slot_layers must be one character or the same length of the assay_layers")
  }

  sc <- import("scanpy", convert = FALSE)
  np <- import("numpy", convert = FALSE)

  obs <- srt@meta.data
  if (ncol(obs) > 0) {
    for (i in seq_len(ncol(obs))) {
      if (is.logical(obs[, i])) {
        obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }

  var <- srt[[assay_X]]@meta.features[features, , drop = FALSE]
  if (ncol(var) > 0) {
    for (i in seq_len(ncol(var))) {
      if (is.logical(var[, i]) && !identical(colnames(var)[i], "highly_variable")) {
        var[, i] <- factor(as.character(var[, i]), levels = c("TRUE", "FALSE"))
      }
    }
  }
  if (length(VariableFeatures(srt, assay = assay_X) > 0)) {
    if ("highly_variable" %in% colnames(var)) {
      var <- var[, colnames(var) != "highly_variable"]
    }
    var[["highly_variable"]] <- features %in% VariableFeatures(srt, assay = assay_X)
  }

  # X <- t(as_matrix(slot(srt@assays[[assay_X]], slot_X)[features, , drop = FALSE]))
  X <- t(GetAssayData(srt, assay = assay_X, slot = slot_X)[features, , drop = FALSE])
  adata <- sc$AnnData(
    X = np_array(X, dtype = np$float32),
    obs = obs,
    var = cbind(data.frame(features = features), var)
  )
  adata$var_names <- features

  layer_list <- list()
  for (assay in names(srt@assays)[names(srt@assays) != assay_X]) {
    if (assay %in% assay_layers) {
      layer <- t(GetAssayData(srt, assay = assay, slot = slot_layers[assay]))
      if (!identical(dim(layer), dim(X))) {
        if (all(colnames(X) %in% colnames(layer))) {
          layer <- layer[, colnames(X)]
        } else {
          stop(
            "The following features in the '", assay_X, "' assay can not be found in the '", assay, "' assay:\n  ",
            paste0(head(colnames(X)[!colnames(X) %in% colnames(layer)], 10), collapse = ","), "..."
          )
        }
      }
      layer_list[[assay]] <- layer
    } else {
      if (isTRUE(verbose)) {
        message("Assay '", assay, "' is in the srt object but not converted.")
      }
    }
  }
  if (length(layer_list) > 0) {
    adata$layers <- layer_list
  }

  reduction_list <- list()
  for (reduction in names(srt@reductions)) {
    reduction_list[[paste0(reduction)]] <- srt[[reduction]]@cell.embeddings
  }
  if (length(reduction_list) > 0) {
    adata$obsm <- reduction_list
  }

  # varm_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   if (ncol(srt[[reduction]]@feature.loadings) > 0) {
  #     varm_list[[paste0(reduction, "_feature.loadings")]] <- srt[[reduction]]@feature.loadings
  #   }
  # }
  # if (length(varm_list) > 0) {
  #   adata$varm <- varm_list
  # }

  obsp_list <- list()
  for (graph in names(srt@graphs)) {
    obsp_list[[graph]] <- srt[[graph]]
  }
  for (neighbor in names(srt@neighbors)) {
    obsp_list[[neighbor]] <- srt[[neighbor]]
  }
  if (length(obsp_list) > 0) {
    adata$obsp <- obsp_list
  }

  uns_list <- list()
  # for (reduction in names(srt@reductions)) {
  #   uns_list[[paste0(reduction,"_stdev")]] <- srt[[reduction]]@stdev
  #   uns_list[[paste0(reduction,"_misc")]] <- srt[[reduction]]@misc
  # }
  # for (image in names(srt@images)) {
  #   uns_list[[image]] <- srt[[image]]
  # }
  if (isTRUE(convert_misc)) {
    for (nm in names(srt@misc)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@misc[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'misc' slot is not converted.")
    }
  }
  if (isTRUE(convert_tools)) {
    for (nm in names(srt@tools)) {
      if (nm != "") {
        uns_list[[nm]] <- srt@tools[[nm]]
      }
    }
  } else {
    if (isTRUE(verbose)) {
      message("'tools' slot is not converted.")
    }
  }
  if (length(uns_list) > 0) {
    adata$uns <- uns_list
  }

  return(adata)
}

#' Convert an anndata object to a seurat object using reticulate
#'
#' @param adata a connected python anndata object.
#'
#' @examples
#' data("pancreas_sub")
#' adata <- srt_to_adata(pancreas_sub)
#' adata <- RunPAGA(adata = adata, group_by = "SubCellType", linear_reduction = "PCA", nonlinear_reduction = "UMAP")
#' srt <- adata_to_srt(adata)
#' srt
#'
#' ### Or convert a h5ad file to Seurat object
#' # library(reticulate)
#' # check_Python("scanpy")
#' # sc <- import("scanpy")
#' # adata <- sc$read_h5ad("pancreas.h5ad")
#' # srt <- adata_to_srt(adata)
#' # srt
#' @importFrom Seurat CreateSeuratObject CreateAssayObject CreateDimReducObject AddMetaData
#' @importFrom SeuratObject as.Graph as.sparse
#' @importFrom reticulate iterate
#' @importFrom Matrix t
#' @export
adata_to_srt <- function(adata) {
  if (!inherits(adata, "python.builtin.object")) {
    stop("'adata' is not a python.builtin.object.")
  }
  sc <- import("scanpy", convert = TRUE)
  np <- import("numpy", convert = TRUE)

  x <- t(py_to_r_auto(adata$X))
  if (!inherits(x, "dgCMatrix")) {
    x <- as.sparse(x[1:nrow(x), , drop = FALSE])
  }
  rownames(x) <- py_to_r_auto(adata$var_names$values)
  colnames(x) <- py_to_r_auto(adata$obs_names$values)

  metadata <- NULL
  if (length(adata$obs_keys()) > 0) {
    metadata <- as.data.frame(py_to_r_auto(adata$obs))
    colnames(metadata) <- make.names(colnames(metadata))
  }

  srt <- CreateSeuratObject(counts = x, meta.data = metadata)

  if (inherits(adata$layers, "python.builtin.object")) {
    keys <- iterate(adata$layers$keys())
  } else {
    keys <- names(adata$layers)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      layer <- py_to_r_auto(adata$layers[[k]])
      if (!inherits(layer, c("Matrix", "matrix"))) {
        stop(paste0("The object in '", k, "' layers is not a matrix: ", paste0(class(adata$layers[[k]]), collapse = ",")))
      }
      layer <- t(layer)
      if (!inherits(layer, "dgCMatrix")) {
        layer <- as.sparse(layer[1:nrow(layer), , drop = FALSE])
      }
      rownames(layer) <- py_to_r_auto(adata$var_names$values)
      colnames(layer) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateAssayObject(counts = layer)
    }
  }

  if (inherits(adata$obsm, "python.builtin.object")) {
    keys <- iterate(adata$obsm$keys())
  } else {
    keys <- names(adata$obsm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsm <- tryCatch(py_to_r_auto(adata$obsm[[k]]), error = identity)
      if (inherits(obsm, "error")) {
        warning("'obsm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsm, "matrix")) {
        obsm <- as_matrix(obsm)
      }
      k <- gsub(pattern = "^X_", replacement = "", x = py_to_r_auto(k))
      colnames(obsm) <- paste0(k, "_", seq_len(ncol(obsm)))
      rownames(obsm) <- py_to_r_auto(adata$obs_names$values)
      srt[[py_to_r_auto(k)]] <- CreateDimReducObject(embeddings = obsm, assay = "RNA", key = paste0(gsub(pattern = "_", replacement = "", x = k), "_"))
    }
  }

  if (inherits(adata$obsp, "python.builtin.object")) {
    keys <- iterate(adata$obsp$keys())
  } else {
    keys <- names(adata$obsp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      obsp <- tryCatch(py_to_r_auto(adata$obsp[[k]]), error = identity)
      if (inherits(obsp, "error")) {
        warning("'obsp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(obsp, "dgCMatrix")) {
        obsp <- as.sparse(obsp[1:nrow(obsp), , drop = FALSE])
      }
      colnames(obsp) <- py_to_r_auto(adata$obs_names$values)
      rownames(obsp) <- py_to_r_auto(adata$obs_names$values)
      obsp <- as.Graph(obsp[seq_len(nrow(obsp)), , drop = FALSE])
      DefaultAssay(object = obsp) <- "RNA"
      srt[[py_to_r_auto(k)]] <- obsp
    }
  }

  if (length(adata$var_keys()) > 0) {
    srt[["RNA"]] <- AddMetaData(srt[["RNA"]], metadata = as.data.frame(py_to_r_auto(adata$var)))
  }

  if (inherits(adata$varm, "python.builtin.object")) {
    keys <- iterate(adata$varm$keys())
  } else {
    keys <- names(adata$varm)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varm <- tryCatch(py_to_r_auto(adata$varm[[k]]), error = identity)
      if (inherits(varm, "error")) {
        warning("'varm: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varm, "matrix")) {
        varm <- as_matrix(varm)
      }
      colnames(varm) <- paste0(py_to_r_auto(k), "_", seq_len(ncol(varm)))
      rownames(varm) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.loadings"]][[py_to_r_auto(k)]] <- varm
    }
  }

  if (inherits(adata$varp, "python.builtin.object")) {
    keys <- iterate(adata$varp$keys())
  } else {
    keys <- names(adata$varp)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      varp <- tryCatch(py_to_r_auto(adata$varp[[k]]), error = identity)
      if (inherits(varp, "error")) {
        warning("'varp: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(varp, "matrix")) {
        varp <- as_matrix(varp)
      }
      colnames(varp) <- py_to_r_auto(adata$var_names$values)
      rownames(varp) <- py_to_r_auto(adata$var_names$values)
      srt[["RNA"]]@misc[["feature.graphs"]][[py_to_r_auto(k)]] <- varp
    }
  }

  if (inherits(adata$uns, "python.builtin.object")) {
    keys <- iterate(adata$uns$keys())
  } else {
    keys <- names(adata$uns)
  }
  if (length(keys) > 0) {
    for (k in keys) {
      uns <- tryCatch(py_to_r_auto(adata$uns[[k]]), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      uns <- tryCatch(check_python_element(uns), error = identity)
      if (inherits(uns, "error")) {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
      if (!inherits(uns, "python.builtin.object")) {
        srt@misc[[py_to_r_auto(k)]] <- uns
      } else {
        warning("'uns: ", k, "' will not be converted. You may need to convert it manually.", immediate. = TRUE)
        next
      }
    }
  }
  return(srt)
}


