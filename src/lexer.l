%{
#include <stdio.h>
#include "parser.h"
%}
%%

('[^']*')+	{
	/* quoted string (one or more) */
	return LABEL;
 }
 /* NOTE: it seems that some (older?) versions of Flex don't recognize the '{-}'
  * (set difference) operator. For now, I don't attempt to support these. */
[[:graph:]]{-}[();,:'\[\]]+	{
	/* printable characters except ();,:'[] */
	return LABEL;
 }
"("	{ return O_PAREN; }
")"	{ return C_PAREN; }
";"	{ return SEMICOLON; }
","	{ return COMMA; }
":"	{ return COLON; }
\[[^]]*]	/* ignore comments */ ;
[\t ]+	/* ignore whitespace */ ;
\n 	;

%%

