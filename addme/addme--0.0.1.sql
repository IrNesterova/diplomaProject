CREATE OR REPLACE FUNCTION
addme(text,integer, integer, text, integer, integer) RETURNS integer AS 'MODULE_PATHNAME','addme'
LANGUAGE C STRICT;


