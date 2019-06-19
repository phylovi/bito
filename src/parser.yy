%require "3.2"
%language "c++"
%{
#include <stdio.h>
#include <string.h>

void yyerror(const char *str)
{
        fprintf(stderr,"error: %s\n",str);
}

int yywrap()
{
        return 1;
}

int main()
{
        yyparse();
}

%}

%token O_PAREN COLON COMMA C_PAREN SEMICOLON
%token LABEL
/* %token <sval> LABEL */

%%

tree: /* empty */ {
  printf("\tEmpty tree\n");
YYACCEPT;
}
| node SEMICOLON {
YYACCEPT;
}
;

node: leaf {
  printf("\tfound a leaf\n");
}
| inner_node {
  printf("\tfound a inner node\n");
}
;

inner_node: O_PAREN nodelist C_PAREN {
  printf("\tGoing in a level\n");
}
;

nodelist: node {
  printf("\tThe last entry\n");
}
| nodelist COMMA node {
  printf("\tChugging along\n");
}
;

leaf: LABEL {
}
;
