my_code_path <- '/data/jianzhongxiang/genome_sequence/'

read_config <- function(){
  # 使用正则表达式筛选文件名中包含"config"的txt文件
  file_names <- list.files(path = './', pattern = "config.*\\.txt$", full.names = TRUE)
  
  # 检查是否有符合条件的文件
  if (length(file_names) == 0) {
    stop("没有找到config文件!!!!!!!!!!!!!!!!!!")
  }
  
  # 读取所有符合条件的txt文件
  for (file in file_names) {
    data <- read.delim(file, sep = "\t", header = TRUE)
    cat("文件名:", basename(file), "\n")
  }
  return(data)
}

sub_rename <- function(config){
  sample_name_row <- which(config[,1] == "[sample_name]")
  sample_group_row <- which(config[,1] == "[sample_group]")
  
  result <- config[(sample_name_row + 1):(sample_group_row - 1), ]
  return(result)
}

sub_group <- function(config){
  sample_group_row <- which(config[,1] == "[sample_group]")
  result <- config[(sample_group_rows[1] + 1):nrow(config), ]
}

merge_q30 <- function(all_file_name) {
  # 初始化一个空列表用于存储每个文件的数据框
  my_q30 <- list()
  
  # 遍历文件名列表，读取每个文件并存储到列表中
  for (i in seq_along(all_file_name)) {
    file_path <- paste0(all_file_name[i], '.q30.txt')  # 构造完整的文件路径
    if (file.exists(file_path)) {  # 检查文件是否存在
      my_q30[[i]] <- read.delim(file = file_path, header = TRUE, check.names = TRUE)
    } else {
      warning(paste0("文件 ", file_path, " 不存在，已跳过。"))
    }
  }
  
  # 使用rbind将所有数据框合并成一张总表
  if (length(my_q30) > 0) {
    total_table <- do.call(rbind, my_q30)
    write.table(total_table,file = paste0('q30.txt'),sep = '\t',row.names = F)
  } else {
    stop("请检查是否有q30文件！！！！！！！！！")
  }
}

gzip_multiplex <- function(folder_path) {
  # 检查文件夹是否存在，若不存在则抛出错误
  if (!dir.exists(folder_path)) {
    stop("Folder does not exist")
  }
  
  # 加载parallel包以支持并行处理
  library(parallel)
  
  # 获取文件夹中所有文件的完整路径
  files <- list.files(folder_path, full.names = TRUE)
  
  # 使用file.info过滤出普通文件（排除目录和无效条目）
  file_info <- file.info(files)
  files <- files[!is.na(file_info$isdir) & !file_info$isdir]
  
  # 并行执行gzip压缩，使用所有可用核心
  mclapply(files, function(file) system(paste("gzip", shQuote(file))), mc.cores = 8)
}

uncompress_multiplex <- function(data_dir = './raw_dir', max_threads = 8, uncompress_dir = './uncompressed') {
  # 检查依赖包是否安装
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' required but not installed")
  }
  
  # 记录开始时间
  message("uncompress_multiplex_cpu start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message("Parameters:\n",
          "Data directory: ", data_dir, "\n",
          "Max threads: ", max_threads, "\n",
          "Uncompress dir: ", uncompress_dir)
  
  # 检查输入目录是否存在
  if (!dir.exists(data_dir)) {
    stop("Data directory does not exist: ", data_dir)
  }
  
  # 创建输出目录
  if (!dir.exists(uncompress_dir)) {
    dir.create(uncompress_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 获取所有 .gz 文件
  gz_files <- list.files(data_dir, pattern = "\\.gz$", full.names = TRUE)
  if (length(gz_files) == 0) {
    stop("No .gz files found in: ", data_dir)
  }
  
  # 打印待处理文件
  message("\nFiles to process (", length(gz_files), "):")
  message(paste(head(basename(gz_files), 5), collapse = "\n"))
  if (length(gz_files) > 5) message("... and ", length(gz_files) - 5, " more")
  
  # 并行处理函数
  process_file <- function(gz_file) {
    file_name <- sub("\\.gz$", "", basename(gz_file))  # 移除 .gz 扩展名
    out_file <- file.path(uncompress_dir, file_name)
    
    if (file.exists(out_file)) {
      message("File exists, skipping: ", out_file)
      return(TRUE)
    }
    
    # 构建解压命令
    cmd <- sprintf("gzip -dc '%s' > '%s'", gz_file, out_file)
    message("Processing: ", basename(gz_file))
    
    # 执行解压命令
    start_time <- Sys.time()
    exit_code <- system(cmd)
    duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
    
    if (exit_code != 0) {
      warning("Failed to decompress: ", basename(gz_file), " (", duration, "s)")
      return(FALSE)
    }
    
    message("Successfully decompressed: ", basename(gz_file), " (", duration, "s)")
    return(TRUE)
  }
  
  # 并行执行
  results <- parallel::mclapply(gz_files, process_file,
                                mc.cores = max_threads,
                                mc.preschedule = FALSE)
  
  # 处理结果统计
  success_count <- sum(unlist(results))
  failure_count <- length(gz_files) - success_count
  
  message("\nProcessing complete: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message("Successfully processed: ", success_count)
  message("Failed files: ", failure_count)
  
  invisible(list(
    total_files = length(gz_files),
    success = success_count,
    failures = failure_count
  ))
}

rename_multiplex <- function(file_path = './uncompressed',config) {
  sample <- sub_rename(config)
  
  files <- list.files(path = './uncompressed', full.names = TRUE)
  
  # 遍历每个文件并重命名
  for (file in files) {
    for (i in 1:nrow(sample)) {
      filename <- basename(file)
      new_filename <- sub(sample[i,1], sample[i,2], filename)
      new_file <- file.path(file_path, new_filename)
      file.rename(file, new_file)
    }
  }
  print('文件重命名已完成！！！！！！！！！！！！！')
}

calculate_q30_multiplex<- function(data_dir = './uncompressed/',phred = 33){
  if(file.exists(file.path(getwd(), "q30.txt"))){
    return('q30.txt已存在！')
  }
  
  get_file_name <- function(path = ".", recursive = FALSE) {
    files <- list.files(path, full.names = FALSE, recursive = recursive)
    files <- files[!file.info(file.path(path, files))$isdir]
    
    # 确保文件名包含扩展名
    files <- files[grepl("\\.[^.]+$", files)]
    
    return(files)
  }
  
  all_file_name <- get_file_name(data_dir)
  
  for (i in all_file_name) {
    print(paste0('开始计算',i,'的q30!'))
    my_cmd <- paste0("python ",my_code_path,"010.calculate_q30.py ", data_dir, i, " --phred ",phred)
    system(my_cmd)
  }
  
  merge_q30(all_file_name)
  print('q30已完成！！！！！！！！！！！！！！！')
}

fastqc_multiplex <- function(data_dir, max_threads=8, contaminants='/data/panzhong/genome_sequence/fastqc/contaminant_list.txt') {
  library(parallel)
  
  # 输出目录逻辑优化
  fastqc_dir <- switch(data_dir,
                       './uncompressed' = './fastqc',
                       './cutadapt' = './fastqc_cutadapt',
                       stop('目录参数异常，请检查脚本运行环境！')
  )
  
  if (!dir.exists(fastqc_dir)) dir.create(fastqc_dir, recursive=TRUE)
  
  input_files <- c(
    list.files(data_dir, pattern="\\.fastq$", full.names=TRUE),
    list.files(data_dir, pattern="\\.fastq\\.gz$", full.names=TRUE)
  )
  
  base_cmd <- paste0("fastqc --quiet --noextract -f fastq",if (!is.null(contaminants)) paste0(" -c ", contaminants)," -o ", fastqc_dir)
  
  # 并行处理逻辑优化
  cl <- makeCluster(min(max_threads, detectCores()))
  clusterExport(cl, "base_cmd", envir=environment())
  
  parLapply(cl, input_files, function(f) {
    output_flag <- file.path(fastqc_dir, paste0(tools::file_path_sans_ext(basename(f)), "_fastqc.zip"))
    if (!file.exists(output_flag)) {
      system2(command=strsplit(base_cmd, " ")[[1]], args=c(shQuote(f), "2>/dev/null"))
    }
  })
  
  stopCluster(cl)
  message("FastQC报告已生成至: ", fastqc_dir)
  print('fastqc已完成！！！！！！！！！！！！！')
}

# 使用示例：run_cutadapt(data_dir = "./uncompressed",cut_dir = "./cutadapt",read_type = "P",max_threads = 4,
#cutadapt_mirna = mirna_samples,cutadapt_ht = ht_samples)
# 定义特殊样本（使用命名列表）
#mirna_samples <- list("sample1" = "CUSTOM_ADAPTOR_FOR_miRNA")
#ht_samples <- list("sample2" = "CUSTOM_ADAPTOR_FOR_HT")
run_cutadapt <- function(data_type = 'cuttag',data_dir = "./uncompressed",cut_dir = "./cutadapt",
                         read_type = "P",
                         file_type = "fastq",
                         adaptor1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",                        
                         adaptor2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
                         phred = 33,
                         max_threads = 8,
                         cutadapt_mirna = NULL,
                         cutadapt_ht = NULL)
{
  library(parallel)
  library(stringr)
  if(data_type == 'cuttag'){
    adaptor1 = adaptor2 = 'CTGTCTCTTATA'
  }
  read_suffix1 = "_R1"
  read_suffix2 = "_R2"
  
  # 创建输出目录
  if (!dir.exists(cut_dir)) dir.create(cut_dir, recursive = TRUE)
  
  # 获取所有fastq文件（包括压缩文件）
  all_files <- list.files(
    path = data_dir,
    pattern = paste0("\\.", file_type, "(\\.gz)?$"),
    full.names = FALSE
  )
  
  # 提取基础样本名
  extract_basename <- function(f) {
    if (read_type == "P") {
      str_remove(f, paste0("(", read_suffix1, "|", read_suffix2, ")\\..+$"))
    } else {
      str_remove(f, paste0("\\.", file_type, "(\\.gz)?$"))
    }
  }
  
  # 获取唯一样本列表
  samples <- unique(sapply(all_files, extract_basename))
  
  # 处理每个样本
  process_sample <- function(sample) {
    tryCatch({
      # 获取完整文件路径
      get_files <- function(pattern) {
        fs <- list.files(
          path = data_dir,
          pattern = paste0("^", sample, pattern, "\\.", file_type, "(\\.gz)?$"),
          full.names = TRUE
        )
        if (length(fs) == 0) return(NA)
        fs
      }
      
      if (read_type == "P") {
        # 处理双端数据
        r1 <- get_files(read_suffix1)
        r2 <- get_files(read_suffix2)
        if (any(is.na(c(r1, r2)))) stop("Missing paired files")
        
        # 检查输出文件
        out_r1 <- file.path(cut_dir, paste0(sample, read_suffix1, ".", file_type))
        out_r2 <- file.path(cut_dir, paste0(sample, read_suffix2, ".", file_type))
        if (file.exists(out_r1) && file.exists(out_r2)) return()
        
        # 构建参数
        base_args <- paste(
          "-m 20 -q 15",
          if (phred == 64) "--quality-base=64",
          "-a", adaptor1,
          "-A", adaptor2
        )
        
        # 特殊处理逻辑
        if (!is.null(cutadapt_mirna) && sample %in% names(cutadapt_mirna)) {
          base_args <- paste(base_args, 
                             "-g CTCGTATGCCGTCTTCTGCTTG",
                             "-G TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAA",
                             "-G AGTTCTGATAACCCACTACCATCGGACCAGCC"
          )
          adaptor2 <- cutadapt_mirna[[sample]]
        } else if (!is.null(cutadapt_ht) && sample %in% names(cutadapt_ht)) {
          base_args <- paste(base_args, "-U 3")
          adaptor2 <- cutadapt_ht[[sample]]
        }
        
        cmd <- sprintf("cutadapt %s -o %s -p %s %s %s",
                       base_args, out_r1, out_r2, r1, r2)
        
      } else {
        # 处理单端数据
        r1 <- get_files("")
        if (is.na(r1)) stop("Missing single-end file")
        
        out_r1 <- file.path(cut_dir, paste0(sample, ".", file_type))
        if (file.exists(out_r1)) return()
        
        base_args <- paste(
          "-m 15 -q 15",
          if (phred == 64) "--quality-base=64",
          "-a", adaptor1
        )
        
        if (!is.null(cutadapt_mirna) && sample %in% names(cutadapt_mirna)) {
          base_args <- paste(base_args, "-g CTCGTATGCCGTCTTCTGCTTG")
        }
        
        cmd <- sprintf("cutadapt %s -o %s %s",
                       base_args, out_r1, r1)
      }
      
      # 执行命令
      message("Processing: ", sample)
      message("[CMD] ", cmd)
      system(cmd)
      
    }, error = function(e) {
      message("Error processing ", sample, ": ", e$message)
    })
  }
  
  # 并行执行
  mclapply(samples, process_sample, mc.cores = max_threads)
  
  message("Cutadapt processing completed!")
}




