# CoNQuiStAdoR
Cell of Nucleic seQuencing uh it's Single to Acid do of Ribo



remember to cite scanpy https://github.com/scverse/scanpy 

**SCANPY: large-scale single-cell gene expression data analysis**

F. Alexander Wolf, Philipp Angerer, Fabian J. Theis

*Genome Biology* 2018 Feb 06. doi: [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0).

# Ideas

## Compression Filter

- indices: bit array of all the non-zero indices
  - int size: # nonzero
- body: bit array of nvz's until there's a gap
  - e.g. body[0] = # of 1s, body[1] = # of 2s, etc
  - int size: # of 1s (a bit wasteful)
- tail: 2 bit arrays, one is nvzs, one is counts
  - e.g. `tail[0]` = [484, 2]
  - tail[1] = [486, 3]
  - etc
  - usually small
  - tail[0]: max nvz
  - tail[1]: max nvz tail count (usually very small)
- total bits
  - B(X): number of bits in X
  - N: number of nonzero values
  - O: number of ones
  - K: size of body
  - T: size of tail
  - M = Max(Nz): maximum nonzero value
  - S: max nzv tail count
  - C: constant extra space (array markers, etc)
  - $$N * B(N) + K * B(O) + T * (B(M) + B(S)) + C$$