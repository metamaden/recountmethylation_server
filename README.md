# recount-methylation-server
Server code for hosting and maintaining Recount-Methylation database.

## Workflow Diagrams
This descriptions and images provide an overview of the various processes handled by Recount Methylation, including processes managed by the server, job queue, and file metadata database, as well as processes that prepare files for analysis.

### Server Process
Recount Methylation includes extensive server functionality that interfaces with the Gene Expression Omnibus to ensure files are properly downloaded, updated, and tracked. Versioning is performed using NTP timestamps. Files are periodically obtained and validated against existing versions of files in the database. Celery server queue management system is used to coordinate the job queue and broker-mediated streaming to a backend SQLite database. Jobs are formed around GSE IDs and filtered GSM lists including valid methylation array idat files. Finally, files are stored in the 'recount-methylation-files' directory, and file metadata is updated in the Recount Methylation MongoDB (RMDB) after a job completes.
![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/server_workflow.tiff "Recount Methylation Server Process")

### Preprocess Workflow
The preprocess workflow utilizes information on the latest versions of files, stored in the Recount Methylation MongoDB (RMDB). RMDB files are validated by searching the root files directory ('./recount-methylation-files'), and metadata and idat files are paired under their corresponding shared GSM ID. Finally, this information is used to form a data sheet including valid file names and paths, and a compressed/'flattened' version of the MetaSRA-pipeline mapped sample metadata. 
![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/preprocess_workflow.tiff "Recount Methylation Server Process")
