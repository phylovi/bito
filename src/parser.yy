/*
Copyright 2019 libsbn project contributors.
libsbn is free software under the GPLv3; see LICENSE file for details.

Based on
https://www.gnu.org/software/bison/manual/html_node/Calc_002b_002b-Parser.html#Calc_002b_002b-Parser
and
https://github.com/tjunier/newick_utils/blob/master/src/newick_parser.y

*** Section: prologue and Bison declarations.
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
  #include "sugar.hpp"
  #include "tree.hpp"
  class Driver;
}

// The parsing context.
%param { Driver& drv }

// Enable location tracking for the parser.
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
  COLON      ":"
  SEMICOLON  ";"
  LPAREN     "("
  RPAREN     ")"
;

%token <std::string> LABEL "label"
%token <std::string> QUOTED "quoted"
%type  <Node::NodePtr> node
%type  <Node::NodePtr> fancy_node
%type  <std::string> leaf
%type  <Node::NodePtr> inner_node
%type  <Node::NodePtrVecPtr> node_list
%type  <std::shared_ptr<Tree>> tree

%printer { yyo << $$; } <*>;

%%
// Grammar rules: how to construct each nonterminal symbol from its parts.

%start tree;
tree:
  fancy_node ";" {
    drv.latest_tree_ = std::make_shared<Tree>($1, drv.branch_lengths_);
    drv.taxa_complete_ = true;
    drv.branch_lengths_.clear();
  };

fancy_node:
  node
| node ":" "label" {
  $$ = $1;

  try {
      SafeInsert(drv.branch_lengths_, $1->Tag(), std::stod($3));
  } catch (...) {
    Failwith("Float conversion failed on branch length '" + $3 +"'");
  }
}

node:
  leaf {
    if (!drv.taxa_complete_) {
      // This is our first tree, so we're going to initialize the taxon set.
      drv.taxa_[$1] = drv.next_id_;
      $$ = Node::Leaf(drv.next_id_);
    drv.next_id_++;
    }
    else {
      // This is not our first tree, so we're going to get taxon numberings from drv.taxa_.
      auto leaf_id = drv.taxa_.find($1);
      if(leaf_id == drv.taxa_.end()) { // leaf $1 not found in taxa_
        Failwith("Taxon '" + $1 + "' did not appear in the first tree.\n" +
         "We only parse lists of trees on the same taxa.");
      }
      $$ = Node::Leaf(leaf_id->second);
    }
  }
| inner_node

leaf:
  "label" {
    $$ = $1;
  }
| "quoted" {
    $$ = $1;
  }

inner_node:
  "(" node_list ")" {
    $$ = Node::Join(std::move(*$2));
  }

node_list:
  fancy_node {
    $$ = std::make_shared<Node::NodePtrVec>();
    $$->push_back($1);
  }
| node_list "," fancy_node {
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
