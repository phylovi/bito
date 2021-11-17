%{
/*
Copyright 2019-2021 bito project contributors.
bito is free software under the GPLv3; see LICENSE file for details.

Based on
https://www.gnu.org/software/bison/manual/html_node/Calc_002b_002b-Scanner.html#Calc_002b_002b-Scanner
and
https://github.com/tjunier/newick_utils/blob/master/src/newick_scanner.l

*** Section: definitions. */
# include <cerrno>
# include <climits>
# include <cstdlib>
# include <cstring> // strerror
# include <string>
# include "driver.hpp"
# include "parser.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"

// Pacify warnings in yy_init_buffer (observed with Flex 2.6.4)
// and GCC 6.4.0, 7.3.0.
#if defined __GNUC__ && !defined __clang__ && 6 <= __GNUC__
# pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

// Of course, when compiling C as C++, expect warnings about NULL.
#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
# pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#define FLEX_VERSION (YY_FLEX_MAJOR_VERSION * 100 + YY_FLEX_MINOR_VERSION)
#if FLEX_VERSION < 206
#error "We require flex version > 2.6."
#endif

/*
After the %option comes the pattern matching.

These are matched greedily (longest match gets priority) and if there are ties
then the first wins.
Note concerning :graph: in code below: this is all non-space printable characters
https://en.cppreference.com/w/cpp/string/byte/isgraph
while :print: is all printable characters
https://en.cppreference.com/w/cpp/string/byte/isprint
also
https://ftp.gnu.org/old-gnu/Manuals/flex-2.5.4/html_mono/flex.html#SEC7
*/
%}

%option noyywrap nounput noinput batch debug

LABEL [[:graph:]]{-}[();,:'\[\]]+
QUOTED ('[^']*')+
BRACKETED_WITH_AMPERSAND \[&[[:print:]]{-}[\[\]]+\]
BLANK [ \t\r]

%{
  // Code run each time a pattern is matched.
  # define YY_USER_ACTION  loc.columns (yyleng);
%}

%%
%{
/* *** Section: rules. */
  // A handy shortcut to the location held by the driver.
  yy::location& loc = drv.location_;
  // Code run each time yylex is called.
  loc.step ();
%}

{BLANK}+   loc.step ();
\n+        loc.lines (yyleng); loc.step ();

"," return yy::parser::make_COMMA(loc);
":" return yy::parser::make_COLON(loc);
";" return yy::parser::make_SEMICOLON(loc);
"(" return yy::parser::make_LPAREN(loc);
")" return yy::parser::make_RPAREN(loc);
{LABEL} return yy::parser::make_LABEL(yytext, loc);
{QUOTED} return yy::parser::make_QUOTED(yytext, loc);
{BRACKETED_WITH_AMPERSAND} return yy::parser::make_BRACKETED_WITH_AMPERSAND(yytext, loc);
.         {
            throw yy::parser::syntax_error
              (loc, "invalid character: " + std::string(yytext));
}
<<EOF>>   return yy::parser::make_END (loc);

%%
/* *** Section: user code. It's just regular C++. */

void
Driver::ScanString(const std::string &str) {
  yy_flex_debug = trace_scanning_;
  yy_scan_string(str.c_str());
}

#pragma GCC diagnostic pop