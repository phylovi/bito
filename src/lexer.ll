%{
#include <stdio.h>
#include "parser.h"
%}
%%

('[^']*')+	{
	/* quoted string (one or more) */
	return yy::parser::token::LABEL;
 }
 /* NOTE: it seems that some (older?) versions of Flex don't recognize the '{-}'
  * (set difference) operator. For now, I don't attempt to support these. */
[[:graph:]]{-}[();,:'\[\]]+	{
	/* printable characters except ();,:'[] */
	return yy::parser::token::LABEL;
 }
"("	{ return yy::parser::token::O_PAREN; }
")"	{ return yy::parser::token::C_PAREN; }
";"	{ return yy::parser::token::SEMICOLON; }
","	{ return yy::parser::token::COMMA; }
":"	{ return yy::parser::token::COLON; }
\[[^]]*]	/* ignore comments */ ;
[\t ]+	/* ignore whitespace */ ;
\n 	;

%%

