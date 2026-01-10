
## Panstarrs
I use queries of the form below except changing the ra and dec coords for all fields.
These queries are submitted to panstarrs dr2 from MAST https://mastweb.stsci.edu/ps1casjobs/SubmitJob.aspx

```
select * into mydb.stds1 from fGetNearbyObjEq(232.664292907031 , 5.98374439399588, 10.0) nb
inner join ObjectThin o on o.objid=nb.objid and o.nDetections>1
inner join MeanObject m on o.objid=m.objid and o.uniquePspsOBid=m.uniquePspsOBid
```
