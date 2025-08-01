# MathX 1.9 Beta - Advanced Mathematical Calculator

## Overview
MathX is a high-performance, feature-rich mathematical calculator designed for students, educators, and professionals. It combines traditional calculation capabilities with advanced features like symbolic computation, high-precision arithmetic, and a built-in learning system.

## Key Features

### Core Calculation Engine
- Supports basic arithmetic operations: `+`, `-`, `*`, `/`, `^` (power)
- Handles complex expressions with parentheses
- Variable assignment and management (`let x = 5`)
- Reverse Polish Notation (RPN) mode support

### Mathematical Functions
- **Trigonometric**: `sin` , `cos`, `tan`, `arcsin`, `arccos`, `arctan`
- **Hyperbolic**: `sinh`, `cosh`, `tanh`
- **Logarithmic/Exponential**: `log`, `ln`, `exp`, `exp2`
- **Special Functions**: `Gamma`, `erf`, `sinc`, `lngamma`
- **Unit Conversion**: `atr` (degrees to radians), `rta` (radians to degrees)
- **Multi-argument Functions**: `min`, `max`, `gcd`, `lcm`, `hypot`

### High-Precision Mode
- Arbitrary-precision integer arithmetic
- Supports addition, subtraction, multiplication, division, and exponentiation
- Enabled via `high precision on` command

### Learning System
- Built-in knowledge base covering:
  - Trigonometry
  - Algebra (polynomials, matrices)
  - Calculus (derivatives, integrals)
- Accessible via `study [topic]` and `search [domain]` commands

### Practical Features
- **History System**:
  - Recall previous results with `!n`
  - Save/load history to files
- **Customization**:
  - Set decimal precision (`setprecision n`)
  - Toggle timing display
  - Multilingual interface (English, Chinese, Russian)
- **Visual Output**:
  - Color-coded messages
  - Clean, formatted display

## Usage Examples

### Basic Calculation
```
>> 2 * sin(pi/4) + 3^2
Result = 10.4142135624
```

### Variable Assignment
```
>> let x 5
>> x^2 + 3*x - 2
Result = 38.0000000000
```

### High-Precision Mode
```
>> high precision on
>> 123456789 * 987654321
Result = 121932631112635269
```

### Learning System
```
>> study derivative
[Displays derivative concepts and rules]
```

## Command Reference

| Command                | Description                          |
|------------------------|--------------------------------------|
| `help`                 | Show all available commands         |
| `let x 5`              | Assign value to variable            |
| `list vars`            | Show all variables                  |
| `clear vars`           | Clear all variables                |
| `setprecision 15`      | Set decimal places (1-80)          |
| `rpn on/off`           | Toggle RPN mode                    |
| `timing on/off`        | Toggle calculation timing display  |
| `language en/zh/ru`    | Change interface language          |
| `!3`                   | Recall 3rd calculation result      |
| `!all`                 | Show complete history              |
| `save history file.txt`| Export history to file             |
| `study [topic]`        | Learn about mathematical concepts  |
| `search [domain]`      | List functions in a domain         |
| `high precision on/off`| Toggle high-precision integer mode |
| `clear`                | Clear screen                       |
| `exit`                 | Quit MathX                         |

## Development Status
Current version: 1.9 Beta  
This is an actively developed project with regular updates and feature additions.

Contribute or report issues on our GitHub repository!
