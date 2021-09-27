# contributing to `bito`

Thank you for your contributions!

bito is written in C++17.
The associated Python module, `vip`, is targeting Python 3.7.


## Overall guidelines

We want the code to be:

1. correct, so we write tests
1. efficient in an algorithmic sense, so we consider algorithms carefully
1. clear to read and understand, so we write code with readers in mind and use code standards
1. fast, so we do profiling to find and eliminate bottlenecks
1. robust, so we use immutable data structures and safe C++ practices
1. simple and beautiful, so we keep the code as minimal and DRY as we can without letting it get convoluted or over-technical

Let's all work together to make the code a pleasure to work with.


### Code structure

* Use the [SOLID principles](https://www.youtube.com/watch?v=Ntraj80qN2k&ab_channel=CppCon) as a guide for good design.
* Factor functionality into functions and methods with logical interfaces.
* Functions should be at most a screenful, with a minimum of arguments.
* Early return is fine at the beginning of a function to avoid wrapping the body of the function in a big if statement. Otherwise strive to have a single return statement in the spirit of [structured programming](https://en.wikipedia.org/wiki/Structured_programming#Early_exit).
* Multiple if statements in a code block which test for the same thing is an antipattern (the ["if" considered harmful talk](https://www.youtube.com/watch?v=z43bmaMwagI&ab_channel=DevWeekEvents) is thought-provoking but a little exaggerated).


### Code clarity

* Use informative and accurate variable and function names.
* Provide sufficient comments and/or documentation to be able to understand how components should be used.
* Prefer [descriptive identifier names and simple coding practices](https://blog.codinghorror.com/coding-without-comments/) to comments documenting implementation.
  If that means having long identifier names, that's fine!
  However, please write documentation if you can't make the code use and operation inherently obvious.
* Prefer clarity over something-that-might-be-slightly-faster (except for the very performance critical core).
* DRY, but prefer clarity over compactness.
* Document methods and functions in header files with comments, unless the method/function name makes its purpose completely obvious.
* If a member variable has an accessor, document the accessor.
* TODO comments don't get merged into master. Rather, make an issue on GitHub and tag future work do with the issue number (e.g. `#262`).


### Use of C++

The [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines) are the authority for how to write C++, and we will follow them.
More generally, we use [clang-tidy](https://clang.llvm.org/extra/clang-tidy/) to check our code according to the `.clang-tidy` file in the root of the repo.
For issues not covered by these guidelines (especially naming conventions), we will use the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
However, we do not use the `k` prefix for consts.
However, the Core Guidelines take priority when these guides differ, such as concerning [passing non-const parameters by reference](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Rf-inout).

We tend towards a modern functional-ish style of C++.

#### Keep inheritance simple

* Limit the depth of inheritance hierarchies.
* If you are considering a pointer to a base class, consider std::variant or [other alternatives](https://www.youtube.com/watch?v=fwXaRH5ffJM&ab_channel=MUCplusplus).

#### Defensive programming

* Use `const` everywhere possible.
* For initialization of `const` variables that requires some computation, [use a lambda](https://herbsutter.com/2013/04/05/complex-initialization-for-a-const-variable/). A switch statement can be useful.
* The default variable initialization should be `const auto`. Range-for loops should loop over `const auto &` by default. ([But don't use auto to store the results of Eigen expressions](https://eigen.tuxfamily.org/dox/TopicPitfalls.html).)
* Use `.at(...)` rather than `[...]` for vectors and maps. This will raise more errors and has fewer [tricky cases](https://stackoverflow.com/a/10821453/467327).
* Add asserts whenever you can't be 100% sure you know what you're getting.
* Always use curly braces for the body of conditionals and loops, even if they are one line.
* Avoid opening namespaces. Never open in header files. Best to open inside a code block such as a function.

#### Keep ownership clear

* Prefer values over pointers. (Raw and shared pointers are like globals!)
* Don't use `new`. Instead, use `make_unique`, or `make_shared` if shared ownership is really necessary.
* Prefer single ownership as per [RAII](https://en.cppreference.com/w/cpp/language/raii).
* Classic/raw pointers are used as non-owning references. Pass smart pointers only when you want to participate in ownership.

#### Use modern C++

* Use modern C++ idioms, including range-for, `auto`, and structured bindings
* [Follow C++ core guidelines for passing parameters](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#fcall-parameter-passing) (F15-F21)
* Prefer returning variables versus non-const arguments to be modified. Because of return value optimization, this doesn't have a performance penalty.
* Use the [rule of zero](https://en.cppreference.com/w/cpp/language/rule_of_three) unless you really need special behavior.
* Use [scoped enums](https://abseil.io/tips/86) rather than C enums.
* [Know](https://www.youtube.com/watch?v=bXkWuUe9V2I&ab_channel=ACCUConference) and [use STL algorithms](https://www.youtube.com/watch?v=W2tWOdzgXHA&ab_channel=25msr) where appropriate.

#### Conventions

* [Follow Abseil](https://abseil.io/tips/88) regarding assignment: use assignment syntax when initializing directly with the intended literal value, and traditional constructor syntax when the initialization is performing logic.
* Don't implement things in non-template header files unless there is a demonstrable performance advantage. One-line accessors are OK.
* Use `[ ]` to index Eigen vectors, and `( )` to index matrices.
* Use lambdas for little internal subroutines. They get CamelCase.


## [Code smells](https://www.youtube.com/watch?v=f_tLQl0wLUM&ab_channel=CppCon)

* Functions with multiple "sections": these should be broken into individual functions
* bool arguments, especially multiple bool arguments. Should we have a function overload or refactor?


### For further work/discussion

* We should have a tool sort the order of things in header files.
* `size_t` vs `Eigen::Index`: what happens if we compare them?
* Use of "Get" for member getters. Note that sometimes we need something, e.g. if we have a `Type` and we want a `GetType()`.
* Standardize list of abbreviations we use, e.g. `idx`. Shall we have an abbreviation for log likelihoods?
* Standardize `fname` vs `path` vs something else.


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
