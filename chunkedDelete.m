function chunkedDelete(targetQuery, batchSize, extraRelvars, autoDetectFK)

error('chunkedDelete is deprecated. Use del() - it now has chunking as an option based on dj.config(''rowDeleteChunk'') ')