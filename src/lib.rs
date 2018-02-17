#![feature(inclusive_range_syntax)]

//! Implementation of [An all-substrings common subsequence algorithm](https://www.sciencedirect.com/science/article/pii/S0166218X07002727)
//!
//! Given two strings s1 and s2, it is possible construct,
//! in O(|s1|\*|s2|) time and O(|s1|+|s2|) space,
//! a structure that can be queried to find the length of
//! all Longest Common Subsequences between s1 and all possible substrings of s2,
//! each query requiring constant time.
//!
//! Some accessor functions are provided to retrieve the matrices and the vectors defined in the paper
//!
//!
//! # Example
//! ```
//! extern crate alcs;
//! use alcs::Alcs;
//!
//! fn main() {
//!    let a = "word";
//!    let b = "hello world";
//!    let va = a.chars().collect::<Vec<char>>();
//!    let vb = b.chars().collect::<Vec<char>>();
//!    let alcs = Alcs::new(&va, &vb);
//!    for i in 0..b.len() {
//!        for (i, j, cij) in alcs.suffix(i) {
//!            println!(r#"LCS between "{}" and "{}" has length {}"#,a,&b[i..j],cij);
//!        }
//!     }
//! }
//! ```
//!Output:
//! ```
//! LCS between "word" and "h" has length 0
//! LCS between "word" and "he" has length 0
//! LCS between "word" and "hel" has length 0
//! ...
//! LCS between "word" and " world" has length 4
//! LCS between "word" and "w" has length 1
//! LCS between "word" and "wo" has length 2
//! LCS between "word" and "wor" has length 3
//! ...
//! LCS between "word" and "d" has length 1
//!```
//!
//! Also, it is defined a trait that allows to fuzzy search a string:
//! ```
//! extern crate alcs;
//! use alcs::FuzzyStrstr;
//!
//! fn main() {
//!    let tsh = 0.7;
//!    let tests = [
//!        ("he.llo.wor.ld.!", "world"),
//!        ("he.llo.word", "world"),
//!        ("hello world", "word"),
//!        ("hello world", "banana"),
//!    ];
//!    for &(h, n) in &tests {
//!        match h.fuzzy_find_str(n, tsh) {
//!            None => {
//!                println!(r#""{}" does not contain "{}""#, h, n);
//!            }
//!            Some((score, sub)) => {
//!                println!(r#""{}" contains "{}" ("{}") with score {}"#, h, n, sub, score);
//!            }
//!        }
//!    }
//! }
//!```
//!Output:
//!```
//! "he.llo.wor.ld.!" contains "world" ("wor.ld") with score 0.8333333
//! "he.llo.word" contains "world" ("word") with score 0.8
//! "hello world" contains "word" ("world") with score 0.8
//! "hello world" does not contain "banana"
//! ```

use std::cmp::{max, min};

/// Constructs the full matrices i_h and i_v
///
/// It requires O(|a|\*|b|) space and O(|a|\*|b|) time
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let (ih,iv) = alcs::compute_mat_ih_iv(&va,&vb);
///    println!("{:?}\n{:?}",ih,iv);
/// ```
///
pub fn compute_mat_ih_iv<T>(a: &[T], b: &[T]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>)
where
    T: Eq,
{
    let na = a.len();
    let nb = b.len();
    let mut ih = vec![vec![0; nb + 1]; na + 1];
    let mut iv = vec![vec![0; nb + 1]; na + 1];

    for j in 0..=nb {
        ih[0][j] = j;
    }

    for l in 0..=na {
        iv[l][0] = 0;
    }

    for l in 1..=na {
        for j in 1..=nb {
            if a[l - 1] != b[j - 1] {
                ih[l][j] = max(iv[l][j - 1], ih[l - 1][j]);
                iv[l][j] = min(iv[l][j - 1], ih[l - 1][j]);
            } else {
                ih[l][j] = iv[l][j - 1];
                iv[l][j] = ih[l - 1][j];
            }
        }
    }
    (ih, iv)
}

/// Constructs the vector I<sub>G</sub>
///
/// It requires O(|a|+|b|) space and O(|a|*|b|) time
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let ig = alcs::compute_vec_ig(&va,&vb);
///    println!("{:?}",ig);
/// ```
pub fn compute_vec_ig<T>(a: &[T], b: &[T]) -> Vec<usize>
where
    T: Eq,
{
    let na = a.len();
    let nb = b.len();
    let mut ih = vec![vec![0; nb + 1], vec![0; nb + 1]];
    let mut iv = vec![vec![0; nb + 1], vec![0; nb + 1]];

    for j in 0..=nb {
        ih[0][j] = j;
    }

    for l in 1..=na {
        iv[1][0] = 0;
        for j in 1..=nb {
            if a[l - 1] != b[j - 1] {
                ih[1][j] = max(iv[1][j - 1], ih[0][j]);
                iv[1][j] = min(iv[1][j - 1], ih[0][j]);
            } else {
                ih[1][j] = iv[1][j - 1];
                iv[1][j] = ih[0][j];
            }
        }
        ih.swap(0, 1);
        iv.swap(0, 1);
    }
    ih.into_iter().next().unwrap()
}

/// Constructs the vectors D<sub>G</sub><sup>0</sup> and V<sub>G</sub> using I<sub>G</sub>
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let ig = alcs::compute_vec_ig(&va,&vb);
///    let (vg,dg) = alcs::compute_vg_dg_from_ig(&va,&vb,&ig);
///    println!("{:?}\n{:?}\n{:?}",ig,vg,dg);
/// ```
pub fn compute_vg_dg_from_ig<T>(
    a: &[T],
    b: &[T],
    ig: &Vec<usize>,
) -> (Vec<Option<usize>>, Vec<Option<usize>>)
where
    T: Eq,
{
    let na = a.len();
    let nb = b.len();
    let mut vg = vec![None; nb + 1];
    let mut dg = vec![Some(0); na + 1];

    let mut i = 1;

    for j in 1..=nb {
        if ig[j] == 0 {
            dg[i] = Some(j);
            i += 1;
        } else {
            vg[ig[j]] = Some(j);
        }
    }

    for l in i..=na {
        dg[l] = None;
    }
    (vg, dg)
}

/// Constructs the vectors D<sub>G</sub><sup>0</sup> and V<sub>G</sub> using the matrix i<sub>h</sub>
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let (ih,iv) = alcs::compute_mat_ih_iv(&va,&vb);
///    let (ig,vg,dg) = alcs::compute_ig_vg_dg_from_ih_mat(&va,&vb,&ih);
///    println!("{:?}\n{:?}",ih,iv);
///    println!("{:?}\n{:?}\n{:?}",ig,vg,dg);
/// ```
pub fn compute_ig_vg_dg_from_ih_mat<T>(
    a: &[T],
    b: &[T],
    ih: &Vec<Vec<usize>>,
) -> (Vec<usize>, Vec<Option<usize>>, Vec<Option<usize>>)
where
    T: Eq,
{
    let ig = ih[ih.len() - 1].clone();
    let (vg, dg) = compute_vg_dg_from_ig(a, b, &ih[ih.len() - 1]);
    (ig, vg, dg)
}

/// Constructs the vectors I<sub>G</sub>, V<sub>G</sub>, D<sub>G</sub><sup>0<sup>
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let (ig,vg,dg) = alcs::alcs(&va,&vb);
///    println!("{:?}\n{:?}\n{:?}",ig,vg,dg);
/// ```
pub fn alcs<T>(a: &[T], b: &[T]) -> (Vec<usize>, Vec<Option<usize>>, Vec<Option<usize>>)
where
    T: Eq,
{
    let ig = compute_vec_ig(a, b);
    let (vg, dg) = compute_vg_dg_from_ig(a, b, &ig);
    (ig, vg, dg)
}

/// Constructs the matrices i<sub>h</sub>, i<sub>v</sub>, and the vectors I<sub>G</sub>, V<sub>G</sub>, D<sub>G</sub><sup>0</sub>
///
/// ```
///    let a = "word";
///    let b = "hello world";
///    let va = a.chars().collect::<Vec<char>>();
///    let vb = b.chars().collect::<Vec<char>>();
///    let (ih,iv,ig,vg,dg) = alcs::alcs_mat(&va,&vb);
///    println!("{:?}\n{:?}\n{:?}\n{:?}\n{:?}",ih,iv,ig,vg,dg);
/// ```
pub fn alcs_mat<T>(
    a: &[T],
    b: &[T],
) -> (
    Vec<Vec<usize>>,
    Vec<Vec<usize>>,
    Vec<usize>,
    Vec<Option<usize>>,
    Vec<Option<usize>>,
)
where
    T: Eq,
{
    let (ih, iv) = compute_mat_ih_iv(a, b);
    let (ig, vg, dg) = compute_ig_vg_dg_from_ih_mat(a, b, &ih);
    (ih, iv, ig, vg, dg)
}

pub struct Alcs {
    ig: Vec<usize>,
}

impl Alcs {
    /// Constructs a structure able to return all the Longest Common Subsequences between a and each substring of b
    ///
    /// It requires O(|a|+|b|) space and O(|a|*|b|) time
    ///
    /// ```
    ///    let a = "word";
    ///    let b = "hello world";
    ///    let va = a.chars().collect::<Vec<char>>();
    ///    let vb = b.chars().collect::<Vec<char>>();
    ///    let alcs = Alcs::new(&va, &vb);
    ///    // this will actually take O(|b|^2)
    ///    for i in 0..b.len() {
    ///        for (i, j, cij) in alcs.suffix(i) {
    ///            println!("LCS between '{}' and '{}' has length {}",a,&b[i..j],cij);
    ///        }
    ///    }
    ///
    /// ```
    pub fn new<T>(a: &[T], b: &[T]) -> Self
    where
        T: Eq,
    {
        Alcs {
            ig: compute_vec_ig(a, b),
        }
    }
    /// Returns an iterator yielding the LCS length for all substrings starting from position `i`
    ///
    pub fn suffix(&self, pos: usize) -> AlcsIterator {
        AlcsIterator::new(&self, pos)
    }
}

pub struct AlcsIterator<'a> {
    alcs: &'a Alcs,
    i: usize,
    j: usize,
    prev: usize,
}

impl<'a> AlcsIterator<'a> {
    fn new(alcs: &'a Alcs, pos: usize) -> Self {
        AlcsIterator {
            alcs,
            i: pos,
            j: pos + 1,
            prev: 0,
        }
    }
}

impl<'a> Iterator for AlcsIterator<'a> {
    type Item = (usize, usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.j >= self.alcs.ig.len() {
            return None;
        }
        let cur = self.prev + (self.alcs.ig[self.j] <= self.i) as usize;
        self.prev = cur;
        self.j += 1;
        Some((self.i, self.j - 1, cur))
    }
}

fn score(b: &str, a: &str, tsh: Option<f32>) -> (f32, usize, usize) {
    let va = a.chars().collect::<Vec<char>>();
    let vb = b.chars().collect::<Vec<char>>();
    let alcs = Alcs::new(&va, &vb);
    let na = a.len();
    let nb = b.len();
    let many = match tsh {
        None => nb,
        Some(tsh) => (na as f32 / tsh) as usize,
    };
    let mut best = (0., 0, 0);
    for i in 0..nb {
        let mut maxrow = (0., 0, 0);
        for (i, j, cij) in alcs.suffix(i).take(many) {
            let cur = cij as f32 / max(j - i, na) as f32;
            if cur > maxrow.0 {
                maxrow = (cur, i, j);
            }
        }
        if maxrow >= best {
            best = maxrow;
        }
    }
    best
}

pub trait FuzzyStrstr<T: AsRef<str>>: AsRef<str> {
    /// Searches s inside self,
    ///
    /// It returns:
    ///
    /// `Some(score,start,end)` if there is a substring s[start..end] achieving a score that is at least the threshold
    ///
    /// `None` otherwise
    ///
    /// For "not too small" values of tsh,
    ///
    /// it requires O(|self|+|s|) space and O(|self|*|s|) time
    fn fuzzy_find_pos(&self, s: T, tsh: f32) -> Option<(f32, usize, usize)> {
        let s = score(self.as_ref(), s.as_ref(), Some(tsh));
        if s.0 > tsh {
            Some(s)
        } else {
            None
        }
    }

    /// Searches s inside self,
    ///
    /// It returns:
    ///
    /// `Some(score,substring)` if there is a substring achieving a score that is at least the threshold
    ///
    /// `None` otherwise
    ///
    /// For "not too small" values of tsh,
    ///
    /// it requires O(|self|+|s|) space and O(|self|*|s|) time
    fn fuzzy_find_str<'a>(&'a self, s: T, tsh: f32) -> Option<(f32, &'a str)> {
        let r = self.fuzzy_find_pos(s, tsh);
        r.map(|(tsh, start, end)| (tsh, &self.as_ref()[start..end]))
    }

    /// Searches s inside self,
    ///
    /// It returns:
    ///
    /// `true` if there is a substring achieving a score that is at least the threshold
    ///
    /// `false` otherwise
    ///
    /// For "not too small" values of tsh,
    ///
    /// it requires O(|self|+|s|) space and O(|self|*|s|) time
    fn fuzzy_contains(&self, s: T, tsh: f32) -> bool {
        self.fuzzy_find_pos(s, tsh).is_some()
    }
}

impl<S, T> FuzzyStrstr<T> for S
where
    S: AsRef<str>,
    T: AsRef<str>,
{
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
