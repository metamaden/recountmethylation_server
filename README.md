# recount-methylation-server
Server code for hosting and maintaining Recount-Methylation database.

## Workflow Diagrams
Descriptions and images summarizing various processes handled by Recount Methylation. These include processes managed by the server, job queue, and file metadata database, as well as processes to prepare files for analysis.

### Server Process
Recount Methylation includes extensive server functionality. It interfaces with the Gene Expression Omnibus, using edirect query utilities and direct ftp calls, to ensure files are properly downloaded, updated, and tracked. Local file versioning is performed using NTP timestamps. Files are periodically obtained and validated against existing versions of files in the database. The Celery server queue management system coordinates execution of the job queue and job status streaming to a backend SQLite database with the RabbitMQ broker. Jobs are formed around GSE IDs and filtered GSM sample lists, and job functions include download and validation of methylation array idat files and GSE soft experiment metadata files. Files are stored in the 'recount-methylation-files' directory, and on job completion, the Recount Methylation MongoDB (RMDB) is updated.
![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/server_workflow.jpg "Recount Methylation Server Process")

### Preprocess Workflow
The preprocess workflow utilizes information on the latest versions of files, stored in the Recount Methylation MongoDB (RMDB), a database of document compilations containing file metadata. RMDB file records are validated by searching the root files directory ('./recount-methylation-files'), and metadata and idat files are paired under their corresponding shared GSM ID. Finally, this information is used to form a data sheet including valid file names and paths, and a compressed/'flattened' version of the MetaSRA-pipeline mapped sample metadata. 
![alt text](https://github.com/metamaden/recount-methylation-server/blob/master/preprocess_workflow.jpg "Recount Methylation Server Process")
