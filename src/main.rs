use anyhow;
use ndarray::Array2;
use std::cmp::Ordering;
use std::fs;

fn naive_table(text: &str) -> Vec<u32> {
    let text = text.as_bytes();
    assert!(text.len() <= u32::MAX as usize);
    let mut table = vec![0u32; text.len()];
    for (i, element) in table.iter_mut().enumerate() {
        *element = i as u32;
    }
    table.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
    table
}

struct SABTree<'a> {
    nblocks: usize,
    block_size: usize,
    t: usize,
    n: usize,
    a: &'a [u32],
    btree: Array2<u32>,
}

#[derive(Debug)]
struct BTreeResult {
    k: usize, // block
    i: usize, // index within block
    r: u32,   // value
}

impl<'a> SABTree<'a> {
    fn new(n: usize, b: usize, a: &'a [u32]) -> Self {
        let block_size = b;
        let nblocks = (n + block_size - 1) / block_size;
        SABTree {
            nblocks,
            block_size,
            t: 0,
            n,
            a: &a,
            btree: Array2::zeros((nblocks, block_size)),
        }
    }

    fn search_fast(&self, x: &str, s: &str) -> BTreeResult {
        let mut k = 0;
        let mut res = BTreeResult {
            k: 0,
            i: 0,
            r: PLACEHOLDER,
        };
        while k < self.nblocks {
            let mut mask: usize = 1 << self.block_size;
            for i in 0..self.block_size {
                println!(
                    "comparing suffix {} at [{},{}] to query {}",
                    &s[self.btree[[k, i]] as usize..],
                    k,
                    i,
                    x
                );
                mask |= match s[self.btree[[k, i]] as usize..].cmp(x) {
                    Ordering::Greater | Ordering::Equal => {
                        println!(">=");
                        1
                    }
                    Ordering::Less => {
                        println!("<");
                        0
                    }
                } << i;
            }
            println!("{:b}", mask);
            let i = mask.trailing_zeros() as usize;
            if i < self.block_size {
                res.k = k;
                res.i = i;
                res.r = self.btree[[k, i]]
            }
            k = self.get_offset(k, i);
            println!("new k = {}", k);
        }
        res
    }

    fn get_offset(&self, k: usize, i: usize) -> usize {
        k * (self.block_size + 1) + i + 1
    }

    fn get_level(&self, idx: usize) -> usize {
        idx / (self.block_size - 1)
    }

    fn get_parent(&self, x: usize) -> (usize, usize) {
        let l: usize = (((x - 1) as f64) / (self.block_size + 1) as f64).floor() as usize;
        let i = (x - 1) - l * (self.block_size + 1);
        (l, i)
    }

    fn build(&mut self, k: usize) {
        if k < self.nblocks {
            for i in 0..self.block_size {
                self.build(self.get_offset(k, i));
                self.btree[[k, i]] = if self.t < self.n {
                    let r = self.a[self.t];
                    self.t += 1;
                    r
                } else {
                    PLACEHOLDER
                };
            }
            let idx = self.get_offset(k, self.block_size);
            self.build(idx);
        }
    }
}

const PLACEHOLDER: u32 = u32::MAX;

fn main() {
    let s = "the_quick_brown_fox_was_quickly_being_quick$";
    println!("{}", s.len());
    let a = naive_table(&s);
    const B: usize = 4;
    let mut sabtree = SABTree::new(a.len(), B, &a);

    sabtree.build(0);
    println!("a = {:?}", a);
    println!("btree =\n{:?}", &sabtree.btree);

    let p = "i}";
    let p0 = "i";

    let lower_bound = sabtree.search_fast(&p0, &s);
    let r2 = lower_bound.r;
    println!("{}", &s[r2 as usize..]);

    let upper_bound = sabtree.search_fast(&p, &s);
    let r = if upper_bound.k != lower_bound.k {
        println!("lower = {:?}, upper = {:?}", lower_bound, upper_bound);
        //upper_bound.r - 1
        //let t = sabtree.get_offset(upper_bound.k, upper_bound.i);
        sabtree.btree[[upper_bound.k, upper_bound.i - 1]]
    } else {
        upper_bound.r
    };

    println!("lb = {:?}, ub = {:?}", lower_bound, upper_bound);

    println!("search for {:?} in {}, [{}, {}]", p, s, r, r2);
    println!("{}", &s[r as usize..]);
    println!("{}", &s[r2 as usize..]);

    /*
    for k in 0..sabtree.nblocks {
        for i in 0..sabtree.block_size {
            let child = sabtree.get_offset(k, i);
            if child < sabtree.nblocks {
                let parent = sabtree.get_parent(child);
                println!(
                    "AT [{}, {}], child is {}, parent is {:?}",
                    k, i, child, parent
                );
            } else {
                println!(
                "AT [{}, {}], calculated child at {}, but that's greater than the number of blocks",
                k, i, child
            );
            }
        }
    }
    */

    let w = [2, 3, 4, 6];
    let z = [2, 3, 5, 6];

    println!("{:?}", &w[1..].cmp(&z[1..]));
}
