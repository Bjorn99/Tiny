# Contributing to Tiny

First off, thank you for considering contributing to Tiny!
## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
  - [Issues](#issues)
  - [Pull Requests](#pull-requests)
- [Development Setup](#development-setup)
- [Project Structure](#project-structure)
- [Coding Guidelines](#coding-guidelines)
- [Testing](#testing)
- [Documentation](#documentation)
- [Commit Messages](#commit-messages)
- [Development Process](#development-process)

## Code of Conduct

This project and everyone participating in it is governed by our [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## Getting Started

### Issues

- **Search existing issues** before creating a new one
- **Use issue templates** when available
- **Provide detailed information** about bugs or feature requests
- **One issue per topic** - don't combine multiple issues into one

#### Bug Reports Should Include:
- Python version
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Any error messages
- Code samples if applicable

#### Feature Requests Should Include:
- Clear description of the feature
- Use cases
- Why it would benefit the project
- Possible implementation ideas (optional)

### Pull Requests

1. **Fork the repository**
2. **Create a new branch** from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bug-fix
   ```
3. **Make your changes**
4. **Write or update tests**
5. **Run the test suite**
6. **Commit your changes** (see [Commit Messages](#commit-messages))
7. **Push to your fork**
8. **Open a Pull Request**

## Development Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/Bjorn99/Tiny.git
   cd Tiny
   ```

2. **Install Poetry**
   ```bash
   curl -sSL https://install.python-poetry.org | python3 -
   ```

3. **Install dependencies**
   ```bash
   poetry install
   ```

4. **Activate virtual environment**
   ```bash
   poetry shell
   ```

5. **Install pre-commit hooks**
   ```bash
   pre-commit install
   ```

## Coding Guidelines

- Follow [PEP 8](https://pep8.org/) style guide
- Use type hints for all functions and methods
- Write docstrings for all public modules, functions, classes, and methods
- Keep functions focused and single-purpose
- Maximum line length is 88 characters (Black default)

### Code Style

We use several tools to maintain code quality:
- **Black** for code formatting
- **isort** for import sorting
- **flake8** for style guide enforcement
- **mypy** for type checking

Run them using:
```bash
poetry run black .
poetry run isort .
poetry run flake8 .
poetry run mypy .
```

## Testing

- Write tests for all new features and bug fixes
- Maintain or improve code coverage
- Use pytest for testing

### Running Tests
```bash
# Run all tests
poetry run pytest

# Run with coverage report
poetry run pytest --cov=tiny tests/

# Run specific test file
poetry run pytest tests/test_sequence.py
```

## Documentation

- Keep README.md updated
- Document all new features
- Update example files if needed
- Add docstrings to all public interfaces
- Update type hints

## Commit Messages

Follow the [Conventional Commits](https://www.conventionalcommits.org/) specification:

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

Types:
- feat: New feature
- fix: Bug fix
- docs: Documentation only changes
- style: Changes that don't affect the meaning of the code
- refactor: Code change that neither fixes a bug nor adds a feature
- test: Adding missing tests or correcting existing tests
- chore: Changes to the build process or auxiliary tools

Example:
```
feat(analysis): add support for RNA sequences

- Added RNA sequence validation
- Updated molecular weight calculations
- Added RNA-specific tests

Closes #123
```

## Development Process

1. **Pick an Issue**
   - Comment on the issue to let others know you're working on it
   - Get assigned to the issue

2. **Create a Branch**
   ```bash
   git checkout -b type/description
   ```

3. **Make Changes**
   - Write code
   - Write tests
   - Update documentation

4. **Test Your Changes**
   ```bash
   poetry run pytest
   poetry run black .
   poetry run isort .
   poetry run flake8 .
   ```

5. **Create Pull Request**
   - Fill out the PR template
   - Link related issues
   - Request review

6. **Address Review Comments**
   - Make requested changes
   - Push updates
   - Respond to comments

7. **Merge**
   - Squash and merge once approved
   - Delete branch

## Questions?

Feel free to create an issue or reach out to the maintainers if you have any questions!
