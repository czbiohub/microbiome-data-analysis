################################################################################

                                GTDB-Tk v0.3.2
________________________________________________________________________________

     GTDB-Tk requires ~25G of external data which needs to be downloaded
     and unarchived. This can be done automatically, or manually:

     1. Run the command download-db.sh to automatically download to:
        /opt/conda/share/gtdbtk-0.3.2/db/

     2. Manually download the latest reference data:
        https://github.com/Ecogenomics/GTDBTk#gtdb-tk-reference-data; Or
        s3://czbiohub-microbiome/ReferenceDBs/GTDB/20191008_R04-RS89/gtdbtk_r89_data.tar.gz

     2b. Set the GTDBTK_DATA_PATH environment variable in the files:
         /opt/conda/etc/conda/activate.d
         /opt/conda/etc/conda/deactivate.d

################################################################################

Alternatively, download the db locally, and use the run_gtdb_classify.sh