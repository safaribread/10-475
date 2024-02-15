# 10-475

      # entry Thurs Feb 8, 2024 
      # Lab meeting minutes: 
        # 1) redo metadata files and check and make sure all of them have (smoker status) 
        # 2) and if needed, convert some data types (e.g. if the smoking level is classified in numbers, convert it) 
        # 3) and give each data line a “source” file to cleanly say which study it came from 
        # 4) And finally combine the metadata files 
              # Try to do it on R, but if not, do it on Excel 
        # 5) Run the chime process for the data sets independently - don't trim them all the same. Trim them as you need to
              # ^if the datasets were all single/paired-end reads and their reads were the same length, you could treat them all the same from the import, but unlikely 
        # 6) Use the Qiime command to merge the data files 
        # 7) DataWrangling section on the proposal: describe how the metadata tables were consolidated 

        #entry Thurs Feb 15, 2024
        #Lab meeting agenda: 
        # 1) Talk about compiled metadata file 
        # 2) Task split
        # 3) Demultiplexing and Quiime
        # 4) Solidify research questions

        # Lab meeting minutes: 
        # 1) Went over compiled metadata files
        # 2) Discussed research question (add more detail)
        # 3) Examined sample sizes (400+ in total; EC 35; QC quite low, might need to be skipped)
        # 4) Suggested doing diversity matrix in R instead of Quiime, to factor in different sample types
        # 5) Outline of aimes for proposal (Due 25th Feb): 1. Metadata wrangling (replace all spaces on the sheet with underscores, simplify the column names); 2. Quiime 2 pipeline (analyze each file separately, merge at .qza, make alpha refraction curve) 3. Alpha and Beta diversity, all of the matrix, done in R. Use alpha and beta diversity to find confounding variables (e.g. locations of samples, location of city change to just country) 4. Do a coremicrobiome analysis 5. Do differential abundance analysis
        # 6) Meeting is optional next week
        # 7) Going through proposal outline
