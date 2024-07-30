# How to Contribute
Thank you for showing interest in contributing to `primerForge`! Here are few guidelines for contributing

## Getting Started
* Make sure you have a GitHub account
* **Add an issue to the github site**. Discuss the problem you are fixing with the maintainer.
  * Has the issue already been discussed on the GitHub issues page?
  * Be detailed and focused

## Making changes
* Fork the `primerForge` repo.
* Create your own new branch and name it with something descriptive and your username, e.g., `dr-joe-wirth-issue-21`.
* Fix the issue on your new branch
  * Make commits of logical and atomic units
  * Add any new unit tests to help with longetivity. Unit tests have a file extension `_test.py`, rely on the `unittest` framework in Python3, and should be placed in the directory `bin/unit_tests`.
    - unit tests can be run using the command `python3 -m unittest discover -s ./primerForge/bin/unit_tests/ -p "*_test.py"`
  * If the unit test(s) show any errors on your local computer, then fix these errors.
* Make a pull request __with the `devel` branch__ when all unit tests show no errors.
