%{
#include <stdio.h>
#include <string.h>

int yylex();
int yyparse();

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

%token O_PAREN COMMA C_PAREN SEMICOLON LABEL

%%

tree: /* empty */ {
  printf("\tEmpty tree\n");
}
| node SEMICOLON
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

leaf: LABEL
;
