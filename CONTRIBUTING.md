# Contributing Guidelines

Thank you for your interest in contributing to `MetaboStatR`! I appreciate your efforts to improve this project. These guidelines will help you understand the contribution process and ensure a smooth collaboration.

## Code of Conduct

To foster an open and welcoming environment, I expect all contributors to adhere to my Code of Conduct. Please ensure your interactions are respectful and constructive.

## How to Contribute

I welcome contributions of all kinds, including bug reports, feature requests, documentation improvements, and code enhancements.

### Reporting Bugs

If you find a bug, please help me by reporting it. Before submitting a new bug report, please check the existing issues to see if the bug has already been reported.

When reporting a bug, provide:
* A clear and concise description of the bug.
* Steps to reproduce the behavior.
* Expected behavior.
* Screenshots or error messages, if applicable.
* Your operating system and R version.

### Suggesting Enhancements

I am always open to new ideas. If you have a suggestion for an enhancement, please open an issue to discuss it. Clearly describe the enhancement and its potential benefits to the project.

### Submitting Pull Requests

For code contributions, please follow these steps:

1.  **Fork the Repository**: Fork the `MetaboStatR` repository to your GitHub account.
2.  **Clone Your Fork**: Clone your forked repository to your local machine.
    ```bash
    git clone [https://github.com/jllcalorio/MetaboStatR.git](https://github.com/jllcalorio/MetaboStatR.git)
    cd MetaboStatR
    ```
3.  **Create a New Branch**: Create a new branch for your changes. Use a descriptive branch name.
    ```bash
    git checkout -b feature/your-feature-name
    # or
    git checkout -b bugfix/your-bug-fix
    ```
4.  **Make Your Changes**: Implement your changes or additions.
    * **Coding Style**: Please adhere to the existing coding style and conventions within the project. Maintain consistency to ensure readability and maintainability.
    * **Tests**: If your changes introduce new functionality or fix a bug, please include appropriate unit tests to ensure stability and correctness.
    * **Documentation**: Update any relevant documentation (e.g., Roxygen comments, README, vignettes) to reflect your changes.
5.  **Commit Your Changes**: Write clear, concise, and descriptive commit messages.
    ```bash
    git commit -m "feat: Add new function for data preprocessing"
    # or
    git commit -m "fix: Resolve issue with data loading"
    ```
6.  **Push to Your Fork**: Push your new branch to your forked repository on GitHub.
    ```bash
    git push origin feature/your-feature-name
    ```
7.  **Open a Pull Request**: Go to the original `MetaboStatR` repository on GitHub and open a new pull request from your branch.
    * **Title**: Provide a clear and concise title for your pull request.
    * **Description**: Explain the purpose of your changes, the problem it solves, or the feature it introduces. Reference any related issues.
    * Ensure all automated tests pass.

## Development Setup

To contribute code, you will need a working R environment. Specific details on setting up the development environment, if any, will be provided within the project's main documentation (e.g., `README.md`).

## Community and Support

If you have questions, need clarification, or want to discuss a potential contribution, please open an issue on the GitHub repository. I will do my best to provide timely support.

## License

By contributing to `MetaboStatR`, you agree that your contributions will be licensed under the project's [License](https://github.com/jllcalorio/MetaboStatR/blob/main/LICENSE).
