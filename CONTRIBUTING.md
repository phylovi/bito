# contributing to `libsbn`

Thank you for your contributions!


## Overall guidelines

We want the code to be:

1. correct, so we write tests
1. efficient in an algorithmic sense, so we consider algorithms carefully
1. clear to read and understand, so we write code with readers in mind and use code standards
1. fast, so we do profiling to find and eliminate bottlenecks
1. robust, so we use immutable data structures and safe C++ practices
1. simple and beautiful, so we keep the code as minimal and DRY as we can without letting it get convoluted or over-technical


## Style

The [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) are the authority for how to write C++, and we will follow them.
More generally, we use [clang-tidy](https://clang.llvm.org/extra/clang-tidy/) to check our code according to the `.clang-tidy` file in the root of the repo.
For issues not covered by these guidelines (especially naming conventions), we will use the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
However, the Core Guidelines take priority when these guides differ, such as concerning [passing non-const parameters by reference](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Rf-inout).

There are certainly violations of these guidelines in the code, so fix them when you see them!


## Specific tips

* Prefer a functional style: returning variables versus modifying them in place. Because of return value optimization, this doesn't have a performance penalty.
* [Follow C++ core guidelines for passing parameters](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#fcall-parameter-passing) (F15-F21)
* [RAII](https://en.cppreference.com/w/cpp/language/raii). No `new`.
* Classic/raw pointers are used as non-owning references. Pass smart pointers only when you want to participate in ownership.
* The default variable initialization should be `const auto`. Range-for loops should loop over `const auto &`. ([But don't use auto to store the results of Eigen expressions](https://eigen.tuxfamily.org/dox/TopicPitfalls.html).)
* Prefer [variable names and simple coding practices](https://blog.codinghorror.com/coding-without-comments/) to code comments.
  If that means having long identifier names, that's fine!
  If you can't make the code use and operation inherently obvious, please write documentation.
* TODO comments don't get merged into master. Rather, make an issue on GitHub.
* Always use curly braces for the body of conditionals and loops, even if they are one line.


## Code smells

* Functions with multiple "sections": these should be broken into individual functions


## Formatting

C++ gets formatted using [clang-format](https://clang.llvm.org/docs/ClangFormat.html), and Python gets formatted using [Black](https://black.readthedocs.io/en/stable/) and [docformatter](https://pypi.org/project/docformatter/).
See the Makefile for the invocations.


## Tests

Add a test for every new feature.

* For C++, we use [doctest](https://github.com/onqtam/doctest).
* For Python, we use [pytest](https://docs.pytest.org/en/latest/).


## [Git flow](https://guides.github.com/introduction/flow/)

* Code changes start by raising an issue proposing the changes, which often leads to a discussion
* Make a branch associated with the issue named with the issue number and a description, such as `4-efficiency-improvements` for a branch associated with issue #4 about efficiency improvements
* If you have another branch to push for the same issue (perhaps a fresh, alternate start), you can just name them consecutively `4-1-blah`, `4-2-etc`, and so on
* Push code to that branch
* Once the code is ready to merge, open a [pull request](https://help.github.com/articles/using-pull-requests/)
* Code review on GitHub
* [Squash and merge](https://help.github.com/en/articles/merging-a-pull-request), [closing the issue via the squash and merge commit message](https://help.github.com/articles/closing-issues-via-commit-messages/)
* Delete branch
