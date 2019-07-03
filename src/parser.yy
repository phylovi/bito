/*
Copyright 2019 Matsen group.
libsbn is free software under the GPLv3; see LICENSE file for details.

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
  #include <cassert>
  #include <string>
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
%type  <Node::NodePtr> tree

%printer { yyo << $$; } <*>;

%%
// Grammar rules: how to construct each nonterminal symbol from its parts.

%start tree;
tree:
  fancy_node ";" {
    auto t = std::make_shared<Tree>($1, drv.branch_lengths_);
    //Tree t($1, drv.branch_lengths_);
    drv.latest_tree_ = $1;
    drv.first_tree_ = false;
    drv.branch_lengths_.clear();
  };

fancy_node:
  node
| node ":" "label" {
  $$ = $1;

  try {
      assert(drv.branch_lengths_.insert(
        std::pair<int64_t, double>($1->Tag(), std::stof($3))).second);
  } catch (...) {
    std::cerr << "Float conversion failed on branch length '" << $3 << "'\n'";
    abort();
  }
}

node:
  leaf {
    if (drv.first_tree_) {
      // This is our first tree, so we're going to initialize the taxon set.
      drv.taxa_[$1] = drv.next_id_;
      $$ = Node::Leaf(drv.next_id_);
    drv.next_id_++;
    }
    else {
      // This is not our first tree, so we're going to get taxon numberings from drv.taxa_.
      auto leaf_id = drv.taxa_.find($1);
      if(leaf_id == drv.taxa_.end()) { // leaf $1 not found in taxa_
        std::cout << "Taxon '" << $1 << "' did not appear in the first tree.\n";
        std::cout << "We only parse lists of trees on the same taxa.\n";
        abort();
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
  // TODO think more about this dereferencing of a shared pointer. I think that this gets return value optimized because
  // $2 goes out of scope, but...
    $$ = Node::Join(*$2);
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
