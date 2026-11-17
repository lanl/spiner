.. _contributing-doc:

.. _Catch2: https://github.com/catchorg/Catch2/blob/devel/docs/tutorial.md

Contributing
=============

If you have any trouble with the project, or are interested in
participating, please contact us by creating an issue on the GitHub
repository, or submit a pull request!

Pull request protocol
----------------------

There is a pull request template that will be auto-populated when you
submit a pull request. A pull request should have a summary of
changes. You should also add tests for bugs fixed or new features you
add.

We have a changelog file, ``CHANGELOG.md``. After creating your pull
request, add the relevant change and a link to the PR in the
changelog.

.. note::

  The changelog is a new requirement, added to make spiner more
  in-line with other codes in the singularity-eos ecosystem.

Before a pull request will be merged, the code should be formatted. We
use clang-format for this, pinned to version 10. You can automatically
trigger ``clang-format`` in two ways: first you can run the script
``format.sh``; second you can type ``make
format_spiner`` after configuring the code with ``clang-format``
discoverable by ``cmake``. The former script takes two CLI arguments
that may be useful, ``CFM``, which can be set to the path for your
clang-format binary, and ``VERBOSE``, which if set to ``1`` adds
useful output. For example:

.. code-block:: bash

    CFM=clang-format-19 VERBOSE=1 ./format.sh

Several sets of tests are triggered on a pull request: a static format
check, a docs build, and unit tests of analytic models and the stellar
collapse model. These are run through GitHub's CPU infrastructure. We
have a second set of tests run on a wider set of architectures that
also access the Sesame library, which we are not able to make public.

The docs are built but not deployed on PRs from forks, and the
internal tests will not be run automatically. So when the code is
ready for merge, you must ask a project maintainer to trigger the
remaining tests for you.

AI-assisted coding
-------------------

``spiner`` requires that if AI was used to assist in code
generation, a disclaimer must be made in a comment in the relevant
file. For example:

..code-block:: c++

  // This file was made in part with generative AI.

Also if agentic AI was used, please have your agent dump a
"proposed plan" markdown file in the ``plan_histories`` folder. This
provides an LLM-readable history of machine-generated changes and
helps disentangle human-made choices from machine-made ones. For
example, if you used codex or claude code, use a workflow like this
one:

1. Ask the agentic framework to propose a plan targeting your problem.
2. Tell it to dump the plan into a new file in ``plan_histories``
3. Iterate until you're happy with the code and submit an MR.
4. After submitting the MR, rename the new file to be prefixed by the MR number and commit it.

If you submit code to ``spiner`` you own that code and you are
responsible for understanding it. If code is submitted that the author
does not understand, the author will be asked to resubmit a changeset
that they understand.

Finally, please be cognizant of reviewer time and effort. Agentic AI
can create changesets much faster than a human can review them. When
possible, please break up large changes and refactors into
human-parse-able chunks.

Expectations for code review
-----------------------------

From the perspective of the contributor
````````````````````````````````````````

Code review is an integral part of the development process
for ``spiner``. You can expect at least one, perhaps many,
core developers to read your code and offer suggestions.
You should treat this much like scientific or academic peer review.
You should listen to suggestions but also feel entitled to push back
if you believe the suggestions or comments are incorrect or
are requesting too much effort.

Reviewers may offer conflicting advice, if this is the case, it's an
opportunity to open a discussion and communally arrive at a good
approach. You should feel empowered to argue for which of the
conflicting solutions you prefer or to suggest a compromise. If you
don't feel strongly, that's fine too, but it's best to say so to keep
the lines of communication open.

Big contributions may be difficult to review in one piece and you may
be requested to split your pull request into two or more separate
contributions. You may also receive many "nitpicky" comments about
code style or structure. These comments help keep a broad codebase,
with many contributors uniform in style and maintainable with
consistent expectations across the code base. While there is no
formal style guide for now, the regular contributors have a sense for
the broad style of the project. You should take these stylistic and
"nitpicky" suggestions seriously, but you should also feel free to
push back.

As with any creative endeavor, we put a lot of ourselves into our
code. It can be painful to receive criticism on your contribution and
easy to take it personally. While you should resist the urge to take
offense, it is also partly code reviewer's responsibility to create a
constructive environment, as discussed below.

Expectations of code reviewers
````````````````````````````````

A good code review builds a contribution up, rather than tearing it
down. Here are a few rules to keep code reviews constructive and
congenial:

* You should take the time needed to review a contribution and offer
  meaningful advice. Unless a contribution is very small, limit
  the times you simply click "approve" with a "looks good to me."

* You should keep your comments constructive. For example, rather than
  saying "this pattern is bad," try saying "at this point, you may
  want to try this other pattern."

* Avoid language that can be misconstrued, even if it's common
  notation in the community. For example, avoid phrases like "code
  smell."

* Explain why you make a suggestion. In addition to saying "try X
  instead of Y" explain why you like pattern X more than pattern Y.

* A contributor may push back on your suggestion. Be open to the
  possibility that you're either asking too much or are incorrect in
  this instance. Code review is an opportunity for everyone to learn.

* Don't just highlight what you don't like. Also highlight the parts
  of the pull request you do like and thank the contributor for their
  effort.

General principle for everyone
```````````````````````````````

It's hard to convey tone in text correspondence. Try to read what
others write favorably and try to write in such a way that your tone
can't be mis-interpreted as malicious.


Notes for Contributors on navigating/developing code features
-------------------------------------------------------------


Some notes on style and code architecture
``````````````````````````````````````````

* A major influence on code style and architecture is the
  `ten rules for developing safety-critical code`_, by Gerard Holzmann.
  Safety critical code is code that exists in a context where failure
  implies serious harm. A flight controller on an airplane or
  spacecraft or the microcontroller in a car are examples of
  safety-critical contexts. ``spiner`` is not safety-critical
  but many of the coding habits advocated for by Holzmann produce
  long-lived, easy to understand, easy to parse, and easy to maintain code.
  And we take many of the rules to heart. Here are a few that are most
  relevant to ``spiner``. They have been adapted slightly to
  our context.

    #. Avoid complex flow constructs such as gotos.

    #. All loops must have fixed bounds. This prevents runaway
       code. (Note this implies that as a general rule, one should use
       ``for`` loops, not ``while`` loops. It also implies one should
       keep recursion to a minimum.)

    #. Heap memory allocation should only be performed at
       initialization. Heap memory de-allocation should only be
       performed at cleanup.

    #. Restrict the length of functions to a single printed page.

    #. Restrict the scope of data to the smallest possible.

    #. Use the preprocessor sparingly.

    #. Limit pointer use to a single dereference. Avoid pointers of
       pointers when possible.

    #. Be compiler warning aware. Try to address compiler warnings as
       they come up.

.. _ten rules for developing safety-critical code: http://web.eecs.umich.edu/~imarkov/10rules.pdf

* ``spiner`` is a modern C++ code
  and both standard template library capabilities and template
  metaprogramming are leveraged frequently. This can sometimes make
  parsing the code difficult. If you see something you don't
  understand, ask. It may be it can be refactored to be more simple or
  better documented.

* As a general rule, to avoid accidental division by zero, use the
  ``robust::ratio(x, y)`` function provided in
  ``ports-of-call`` instead of writing ``x /  y``.

Performance portability concerns
`````````````````````````````````

``spiner`` is performance portable, meaning it is designed to
run not only on CPUs, but GPUs from a variety of manufacturers,
powered by a variety of device-side development tools such as Cuda,
OpenMP, and OpenACC. This implies several constraints on code
style. Here we briefly discuss a few things one should be aware of.

* **``ports-of-call`` and portability decorators:** Functions that
  should be run on device needs to be decorated with one of the
  following macros: ``PORTABLE_FUNCTION``,
  ``PORTABLE_INLINE_FUNCTION``,
  ``PORTABLE_FORCEINLINE_FUNCTION``. These macros are imported from
  the `ports-of-call`_ library and resolve to the appropriate
  decorations for a given device-side backend such as Cuda so the code
  compiles correctly. Code that doesn't need to run on device,
  such as class constructors, does not need these decorations.

* **Relocatable device code:** It is common in C++ to split code
  between a header file and an implementation file. Functionality that
  is to be called from within loops run on device should not be split
  in this way. Not all accelerator languages support this and the ones
  that do take a performance hit. Instead implement that functionality
  only in a header file and decorate it with
  ``PORTABLE_INLINE_FUNCTION`` or ``PORTABLE_FORCEINLINE_FUNCTION``.

* **Host and device pointers:** Usually accelerators have different
  memory spaces than the CPU they are attached to. So you need to be
  aware that data needs to be copied to an accelerator device to be
  used. If it is not properly copied, the code will likely crash with
  a segfault. In general scalar data such as a single variable (e.g.,
  ``int x``) can be easily and automatically copied to device and you
  don't need to worry about managing it. Arrays and pointers, however,
  are a different story. If you create an array or point to some
  memory on CPU, then you are pointing to a location in memory on your
  CPU. If you try to access it from your accelerator, your code will
  not behave properly. You need to manually copy data from host to
  device in this case.

* **Shallow copies:** As a general rule, large amount of data stored
  within a ``Databox`` object should have "reference-semantics." This
  means that if you copy a ``Databox`` object, it should always be a
  shallow copy, not a deep copy, unless a deep copy is explicitly
  requested. This is for performance reasons and also to simplify the
  management of data on device.

* **Real:** The ``Real`` datatype is either a single precision or
  double precision floating point number, depending on how
  `ports-of-call`_ is configured. For most floating point numbers use
  the ``Real`` type. However, be conscious that sometimes you will
  specifically need a single or double precision number, in which case
  you should specify the type as built into the language.

.. _ports-of-call: https://lanl.github.io/ports-of-call/main/index.html


How to Make a Release
----------------------

``spiner`` uses *semantic versioning*. A version is written
as ``v[major version].[minor version].[patch number]``. To make a new
release, first make a new pull request where you (1) change the
version number in the ``project`` field of the of the top-level
``CmakeLists.txt`` file and (2) add a new release field to the
``CHANGELOG.md``, moving all the changes listed under ``Current Main``
to that release. Then add empty categories for ``Current
Main``. Typically the branch for this merge request should be called
``v[release number]-rc`` for "release candidate." Make sure that the
full test suite passes for this PR.

After that pull request is merged, go to the ``releases`` tab on the
right sidebar on GitHub, and draft a new release. Set the tag to
``v[release number]``, fill the comment with the changes in the
changelog since the last release, and make the release.

Finally, the Spackages must be updated. To do so, you will need to
submit an MR to the spack `spack-packages`_ repository with an updated
``package.py`` with the new release. LANL developers must also submit
to the ``xcap-spackages`` library.

.. _spack-packages: https://github.com/spack/spack-packages

Continuous Integration
----------------------

``spiiner`` has two continuous integration (CI) systems. A public
facing one via GitHub actions and an LANL internal one through a
GitLab instance. The internal one uses `kessel`_. See the Kessel
documentation for more details.

.. _kessel: https://github.com/lanl/kessel
