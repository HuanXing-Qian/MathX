The current MathX is version 1.8 Beta, and to be honest, it will probably stay in Beta for a long time because I have so many ideas—new ones keep popping up almost every moment that I want to add.

# MathX 1.8 Beta - Advanced Mathematical Calculator

## Overview
MathX is a feature-rich, high-performance mathematical calculator and computing environment developed in C++. Designed for both educational and professional use, it combines traditional calculator functionality with advanced features like arbitrary-precision arithmetic, symbolic computation, and an integrated knowledge system.

## Key Features

### Core Calculation Engine
- **Dual-Precision Mode**: Switch between standard floating-point and high-precision integer arithmetic
- **Comprehensive Function Support**:
  - 30+ mathematical functions (trigonometric, logarithmic, statistical, etc.)
  - Special functions (Gamma, erf, sinc, etc.)
  - Unit conversion (degrees/radians)
- **Constants Library**: Predefined mathematical constants (π, e, γ, φ, √2) with 100+ digit precision
- **Variables System**: User-defined variables with persistence

### Advanced Capabilities
- **Parallel FFT Implementation**: Optimized large-number multiplication using parallelized Fast Fourier Transform
- **Dynamic Precision Control**: Adjustable from 1 to 80 decimal places
- **RPN Mode**: View Reverse Polish Notation translation of expressions
- **History System**: Complete calculation history with save/load functionality

### Unique Functionalities
- **Integrated Knowledge Base**:
  - Multilingual support (English, Chinese, Russian)
  - Contextual help system with fuzzy matching
  - Domain-specific function documentation
- **Performance Metrics**: Execution timing for optimization
- **Developer-Friendly**:
  - ANSI color-coded output
  - GCC-optimized compilation
  - Clean, modular architecture

## Technical Specifications
- **Language**: C++17
- **Dependencies**: Standard Library, Windows.h (for console operations)
- **Optimizations**: 
  - `#pragma GCC optimize("O3,unroll-loops")`
  - SIMD-parallelized FFT
  - Move semantics for large-number operations

## Usage Examples
```plaintext
Basic: 2 * sin(pi/4)
Variables: let x 42 then x^2
Precision: setprecision 20
High-precision: 12345678901234567890 * 9876543210987654321
Knowledge: study derivative
```

## Development Roadmap
- [ ] Symbolic computation engine
- [ ] Matrix operations
- [ ] Graphical plotting
- [ ] Plugin system for custom functions
- [ ] Cross-platform compatibility

## Building
```bash
g++ MathX.cpp -o MathX -O3 -march=native -funroll-loops
```

## Known Limitations
- High-precision mode currently limited to integer operations
- Windows console dependency for color output
- Beta status indicates ongoing API changes

## Contribution Guidelines
As the sole developer currently, I welcome:
- Bug reports with reproducible cases
- Performance optimization suggestions
- Mathematical function implementations
- Documentation improvements

This project reflects my passion for both mathematics and systems programming - where elegant algorithms meet optimized implementation. The "beta" tag represents not instability, but rather the endless stream of improvements I continue to implement.

**Developer**: ShiYuze  
**Version**: 1.8 Beta  
**License**: Proprietary (for now)

By the way, since I prefer putting all functions into a single file, if you only need one language, just download the main MathX file and the language pack (the .txt one!), then place them in the same folder—it'll run directly!
