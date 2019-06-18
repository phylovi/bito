%{
#include <stdio.h>
#include "parser.h"
%}
%%

"("	{ return O_PAREN; }
")"	{ return C_PAREN; }
";"	{ return SEMICOLON; }
","	{ return COMMA; }
('[^']*')+	{ return LABEL; }
[[:graph:]]{-}[();,:'\[\]]+	{ return LABEL; }
\n                      /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;

%%

