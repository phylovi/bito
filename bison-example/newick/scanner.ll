%{ /* -*- C++ -*- */
/* *** Section: definitions. */
# include <cerrno>
# include <climits>
# include <cstdlib>
# include <cstring> // strerror
# include <string>
# include "driver.hh"
# include "parser.hh"

// Pacify warnings in yy_init_buffer (observed with Flex 2.6.4)
// and GCC 6.4.0, 7.3.0.
#if defined __GNUC__ && !defined __clang__ && 6 <= __GNUC__
# pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

// Of course, when compiling C as C++, expect warnings about NULL.
#if defined __clang__
# pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
# pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#define FLEX_VERSION (YY_FLEX_MAJOR_VERSION * 100 + YY_FLEX_MINOR_VERSION)
#if FLEX_VERSION < 206
#error "We require flex version > 2.6."
#endif
%}

%option noyywrap nounput noinput batch debug

%{
  // A number symbol corresponding to the value in S.
  yy::parser::symbol_type
  make_NUMBER (const std::string &s, const yy::parser::location_type& loc);
%}

TAXON         [[:graph:]]{-}[();,:'\[\]]+
QUOTED_TAXON  ('[^']*')+
INT           [0-9]+
BLANK         [ \t\r]

%{
  // Code run each time a pattern is matched.
  # define YY_USER_ACTION  loc.columns (yyleng);
%}

%%
%{
/* *** Section: rules. */
  // A handy shortcut to the location held by the driver.
  yy::location& loc = drv.location;
  // Code run each time yylex is called.
  loc.step ();
%}
{BLANK}+   loc.step ();
\n+        loc.lines (yyleng); loc.step ();

","                  return yy::parser::make_COMMA     (loc);
";"                  return yy::parser::make_SEMICOLON (loc);
"("                  return yy::parser::make_LPAREN    (loc);
")"                  return yy::parser::make_RPAREN    (loc);

{INT}                return make_NUMBER (yytext, loc);
{TAXON}              return yy::parser::make_TAXON (yytext, loc);
{QUOTED_TAXON}       return yy::parser::make_QUOTED_TAXON (yytext, loc);
.                    {
                       throw yy::parser::syntax_error
                         (loc, "invalid character: " + std::string(yytext));
}
<<EOF>>              return yy::parser::make_END (loc);

%%
/* *** Section: user code. It's just regular C++. */

yy::parser::symbol_type
make_NUMBER (const std::string &s, const yy::parser::location_type& loc)
{
  errno = 0;
  long n = std::strtol (s.c_str(), NULL, 10);
  if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE))
    throw yy::parser::syntax_error (loc, "integer is out of range: " + s);
  return yy::parser::make_NUMBER ((int) n, loc);
}

void
driver::scan_begin ()
{
  yy_flex_debug = trace_scanning;
  if (file.empty () || file == "-")
    yyin = stdin;
  else if (!(yyin = fopen (file.c_str (), "r")))
    {
      std::cerr << "cannot open " << file << ": " << strerror(errno) << '\n';
      exit (EXIT_FAILURE);
    }
}

void
driver::scan_end ()
{
  fclose (yyin);
}
