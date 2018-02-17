# Alcs
All Longest Common Substrings and Fuzzy String Search for Rust
 Link to [crates.io](https://crates.io/crates/alcs)


 Implementation of [An all-substrings common subsequence algorithm](https://www.sciencedirect.com/science/article/pii/S0166218X07002727)

 Given two strings s1 and s2, it is possible construct,
 in O(|s1|\*|s2|) time and O(|s1|+|s2|) space,
 a structure that can be queried to find the length of
 all Longest Common Subsequences between s1 and all possible substrings of s2,
 each query requiring constant time.

 Some accessor functions are provided to retrieve the matrices and the vectors defined in the paper


 # All LCS
 ```rust
 extern crate alcs;
 use alcs::Alcs;

 fn main() {
    let a = "word";
    let b = "hello world";
    let va = a.chars().collect::<Vec<char>>();
    let vb = b.chars().collect::<Vec<char>>();
    let alcs = Alcs::new(&va, &vb);
    for i in 0..b.len() {
        for (i, j, cij) in alcs.suffix(i) {
            println!(r#"LCS between "{}" and "{}" has length {}"#,a,&b[i..j],cij);
        }
     }
 }
 ```
Output:
 ```
 LCS between "word" and "h" has length 0
 LCS between "word" and "he" has length 0
 LCS between "word" and "hel" has length 0
 ...
 LCS between "word" and " world" has length 4
 LCS between "word" and "w" has length 1
 LCS between "word" and "wo" has length 2
 LCS between "word" and "wor" has length 3
 ...
 LCS between "word" and "d" has length 1
```
# Fuzzy substring
 Also, it is defined a trait that allows to fuzzy search a string:
 ```rust
 extern crate alcs;
 use alcs::FuzzyStrstr;

 fn main() {
    let tsh = 0.7;
    let tests = [
        ("he.llo.wor.ld.!", "world"),
        ("he.llo.word", "world"),
        ("hello world", "word"),
        ("hello world", "banana"),
    ];
    for &(h, n) in &tests {
        match h.fuzzy_find_str(n, tsh) {
            None => {
                println!(r#""{}" does not contain "{}""#, h, n);
            }
            Some((score, sub)) => {
                println!(r#""{}" contains "{}" ("{}") with score {}"#, h, n, sub, score);
            }
        }
    }
 }
```
Output:
```
 "he.llo.wor.ld.!" contains "world" ("wor.ld") with score 0.8333333
 "he.llo.word" contains "world" ("word") with score 0.8
 "hello world" contains "word" ("world") with score 0.8
 "hello world" does not contain "banana"
 ```
