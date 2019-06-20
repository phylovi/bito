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
  // This code gets inserted into the parser header file.
  #include <string>
  #include "sbn.hpp"
  class Driver;
}

// The parsing context.
%param { Driver& drv }

%locations

%define parse.trace
%define parse.error verbose

%code {
#include "driver.hpp"
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
%token <std::string> QUOTED_TAXON "quoted_taxon"
%token <int> NUMBER "number"
%type  <Node::NodePtr> node
%type  <std::string> leaf
%type  <Node::NodePtr> inner_node
%type  <Node::NodePtrVecPtr> node_list
%type  <Node::NodePtr> tree

%printer { yyo << $$; } <*>;

%%
// Grammar rules: how to construct each nonterminal symbol from its parts.

%start tree;
tree:
  node ";" {
    drv.latest_tree_ = $1;
    drv.first_tree_ = false;
  };

node:
  leaf {
    int leaf_id;
    if (drv.first_tree_) {
      // This is our first tree, so we're going to initialize the taxon set.
      drv.taxa[$1] = drv.next_id_;
      $$ = Node::Leaf(drv.next_id_);
    drv.next_id_++;
    }
    else {
      // This is not our first tree, so we're going to get taxon numberings from drv.taxa.
      auto leaf_id = drv.taxa.find($1);
      if(leaf_id == drv.taxa.end()) { // leaf $1 not found in taxa
        std::cout << "Taxon '" << $1 << "' did not appear in the first tree.\n";
        std::cout << "We only parse lists of trees on the same taxa.\n";
        abort();
      }
      $$ = Node::Leaf(leaf_id->second);
    }
  }
| inner_node

leaf:
  "taxon" {
    $$ = $1;
  }
| "quoted_taxon" {
    $$ = $1;
  }

inner_node:
  "(" node_list ")" {
  // TODO think more about this dereferencing of a shared pointer.
    $$ = Node::Join(*$2);
  }

node_list:
  node {
    $$ = std::make_shared<Node::NodePtrVec>();
    $$->push_back($1);
  }
| node_list "," node {
    $1->push_back($3);
    $$ = $1;
  }

%%
// Epilogue: arbitrary C++.

void
yy::parser::error (const location_type& l, const std::string& m)
{
  std::cerr << l << ": " << m << '\n';
}
