INSERT INTO Research ( InternalPatientID, ResearchProject )
SELECT DISTINCT Patients.InternalPatientID, "1923" AS ResearchProject
FROM BulkUpload20171003 LEFT JOIN Patients ON BulkUpload20171003.PatientID = Patients.PatientID;
