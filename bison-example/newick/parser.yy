/* *** Section: prologue and Bison declarations.
  The prologue is broken up into %code blocks, with an optional qualifier
  to describe where it should go in the resulting source file. */

%skeleton "lalr1.cc" /* -*- C++ -*- */
%require "3.4"
%defines

%define api.token.constructor
%define api.value.type variant
%define parse.assert

%code requires {
  # include <string>
  class driver;
}

// The parsing context.
%param { driver& drv }

%locations

%define parse.trace
%define parse.error verbose

%code {
# include "driver.hh"
}

%define api.token.prefix {TOK_}

// Bison declarations: the names, types, and precedence, of symbols.

%token
  END  0  "end of file"
  COMMA      ","
  SEMICOLON  ";"
  LPAREN     "("
  RPAREN     ")"
;

%token <std::string> TAXON "taxon"
%token <int> NUMBER "number"
%type  <int> node
%type  <int> inner_node
%type  <int> node_list

%printer { yyo << $$; } <*>;

%%
// Grammar rules: how to construct each nonterminal symbol from its parts.

%start tree;
tree:
  node ";" {
    drv.result = $1;
    drv.id_counter = 0; // Reset id_counter to zero.
  };

node:
  "taxon" {
    $$ = 1;
    drv.taxa[$1] = drv.id_counter;
    drv.id_counter++;

  }
| inner_node

inner_node:
  "(" node_list ")" {
    $$ = $2;
    drv.id_counter++;
  }

node_list:
  node
  | node_list "," node {
    $$ = $1 + $3;
  }

%%
// Epilogue: arbitrary C++.

void
yy::parser::error (const location_type& l, const std::string& m)
{
  std::cerr << l << ": " << m << '\n';
}
