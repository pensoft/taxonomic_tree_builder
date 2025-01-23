CREATE UNIQUE INDEX taxonid_idx ON taxon (taxonid);

SELECT taxon.taxonid, jsonb_agg(concat(
	'{"taxonid":"', jtaxon.taxonid, '"', ',"scientificname":"', jtaxon.scientificname, '",','"taxonrank":"', jtaxon.taxonrank, '"}')::json order by jtaxon.id ASC) 
FROM taxon
LEFT JOIN taxon as jtaxon on jtaxon.taxonid = any(taxon.parents) 
where taxon.parentnameusageid is not null --and taxon.taxonid = '622BF'
GROUP BY taxon.taxonid
limit 10000