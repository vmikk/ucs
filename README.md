# `ucs` - UC file proccessing tool

The `ucs` tool reads USEARCH cluster format (UC) files, 
which are tab-separated text files commonly used in clustering and database searches. 
Each line in a UC file represents a record corresponding to an input sequence, 
describing cluster-membership information, alignments, and sequence identities. 
For detailed information on the UC format, refer to the [USEARCH manual](https://drive5.com/usearch/manual/opt_uc.html).

## Usage

Check clustering summary 
(estimates the number of unique query and target sequences, 
as well as the number of queries with identical names mapped to multiple targets):

```bash
ucs -i test.uc.gz -s
```
