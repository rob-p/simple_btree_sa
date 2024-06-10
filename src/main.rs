use ndarray::Array2;
use std::cmp::Ordering;

/// Naive implementation of suffix array construction just for
/// testing purposes.  Taken directly from the
/// implementation in the `suffix` crate.
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

/// Represents a "SuffixArrayB-Tree" (SABTree)
struct SABTree {
    nblocks: usize,
    block_size: usize,
    t: usize,
    n: usize,
    btree: Array2<u32>,
}

#[derive(Debug, Clone)]
struct BTreePointer {
    k: usize, // block
    i: usize, // index within block
    r: u32,   // value
}

impl BTreePointer {
    pub fn advance(&mut self, stack: &mut Vec<BTreePointer>, t: &SABTree) -> bool {
        loop {
            // if we are at the end of the current block, then
            // go to the parent
            if self.i >= t.block_size || t.btree[[self.k, self.i]] == PLACEHOLDER {
                // do we have a parent on the stack? If so, visit it
                if let Some(x) = stack.pop() {
                    self.k = x.k;
                    self.i = x.i;
                    self.r = t.btree[[self.k, self.i]];
                    return true;
                } else {
                    loop {
                        // if we are not at the top level, we can
                        // get the parent, and then visit it. In this
                        // case we are visiting the parent after seeing
                        // its children and so we want to place the
                        // *next* node (parent) at this level on the stack.

                        // try and get the parent
                        if let Some((plevel, pidx)) = t.get_parent(self.k) {
                            self.k = plevel;
                            self.i = pidx;
                            // if the successor of the parent is valid, add that
                            // to the stack
                            if pidx + 1 < t.block_size {
                                stack.push(BTreePointer {
                                    k: self.k,
                                    i: self.i + 1,
                                    r: t.btree[[self.k, self.i + 1]],
                                });
                            }
                            // if the parent itself is valid, we can fill
                            // in the value and return it. If the parent
                            // isn't valid (the offset is too great) we
                            // will loop again and look for it's parent
                            if self.i < t.block_size && t.btree[[self.k, self.i]] != PLACEHOLDER {
                                self.r = t.btree[[self.k, self.i]];
                                break;
                            }
                        } else {
                            // there is no parent; we're done
                            return false;
                        }
                    }
                    return true;
                }
            }

            // if we are not at the end of a block check if the
            // current node has a child. If so, move to that
            let possible_child_level = t.get_offset(self.k, self.i + 1);
            if possible_child_level < t.nblocks {
                self.k = possible_child_level;
                self.i = 0;
                self.r = t.btree[[self.k, self.i]];
                return true;
            } else {
                // the current node has no child so we are in the "simple"
                // case of just advancing i
                self.i += 1;
                if self.i >= t.block_size || t.btree[[self.k, self.i]] == PLACEHOLDER {
                    continue;
                }
                self.r = t.btree[[self.k, self.i]];
                return true;
            }
        }
    }
}

impl PartialEq for BTreePointer {
    /// Test if two BTreePointers point to the same node
    /// (tests only the level and offset, ignoring the value)
    fn eq(&self, other: &Self) -> bool {
        self.k == other.k && self.i == other.i
    }
}

struct SABTreeIterator<'a> {
    t: &'a SABTree,
    end: BTreePointer,
    curr: BTreePointer,
    valid: bool,
    stack: Vec<BTreePointer>,
}

impl<'a> SABTreeIterator<'a> {
    fn new(lower_bound: &BTreePointer, upper_bound: &BTreePointer, tree: &'a SABTree) -> Self {
        Self {
            t: tree,
            end: upper_bound.clone(),
            curr: lower_bound.clone(),
            valid: true,
            stack: Vec::with_capacity(8),
        }
    }
}

impl<'a> Iterator for SABTreeIterator<'a> {
    type Item = BTreePointer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr != self.end && self.valid {
            // yield the current value and then
            // advance the pointer.
            let r = self.curr.clone();
            self.valid = self.curr.advance(&mut self.stack, self.t);
            Some(r)
        } else {
            // There is nothing more to iterate
            None
        }
    }
}

/// The implementations below (particularly build and search) are directly based off of the
/// "Implicit Static B-trees" (i.e. S-trees) examples given on the
/// [Algorithmica website](https://algorithmica.org/en/b-tree).
impl SABTree {
    fn new(n: usize, b: usize) -> Self {
        let block_size = b;
        let nblocks = (n + block_size - 1) / block_size;
        SABTree {
            nblocks,
            block_size,
            t: 0,
            n,
            btree: Array2::zeros((nblocks, block_size)),
        }
    }

    /// Returns true if p points within a leaf node and
    /// false otherwise
    #[allow(dead_code)]
    fn is_leaf(&self, p: &BTreePointer) -> bool {
        self.get_offset(p.k, p.i) >= self.nblocks
    }

    fn search(&self, x: &str, s: &str) -> BTreePointer {
        let mut k = 0;
        let mut res = BTreePointer {
            k: 0,
            i: 0,
            r: PLACEHOLDER,
        };
        while k < self.nblocks {
            let mut mask: usize = 1 << self.block_size;
            for i in 0..self.block_size {
                // NOTE: Arghhh... this is frustrating. We'd like to do our best
                // to avoid a conditional branch here. The "PLACEHOLDER" trick does
                // this when we are comparing integers directly, but when we are
                // comparing suffixes like this, we have to avoid indexing the
                // underlying string out of bounds. Think about how we might
                // avoid a conditional here ... is there a good "placeholder" in the
                // suffix array context?
                let idx = self.btree[[k, i]] as usize;
                mask |= if idx != PLACEHOLDER as usize {
                    match s[idx..].cmp(x) {
                        Ordering::Greater | Ordering::Equal => {
                            println!(">=");
                            1
                        }
                        Ordering::Less => {
                            println!("<");
                            0
                        }
                    }
                } else {
                    1
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

    /// If the current node has a parent, returns Some((parent_level, parent_offset)), otherwise
    /// if the current node is at the root, it returns None.
    fn get_parent(&self, x: usize) -> Option<(usize, usize)> {
        if x >= 1 {
            let l: usize = (((x - 1) as f64) / (self.block_size + 1) as f64).floor() as usize;
            let i = (x - 1) - l * (self.block_size + 1);
            Some((l, i))
        } else {
            None
        }
    }

    fn build(&mut self, a: &[u32], k: usize) {
        if k < self.nblocks {
            for i in 0..self.block_size {
                self.build(a, self.get_offset(k, i));
                self.btree[[k, i]] = if self.t < self.n {
                    let r = a[self.t];
                    self.t += 1;
                    r
                } else {
                    PLACEHOLDER
                };
            }
            let idx = self.get_offset(k, self.block_size);
            self.build(a, idx);
        }
    }
}

const PLACEHOLDER: u32 = u32::MAX;

fn main() -> anyhow::Result<()> {
    let s = "the_quick_brown_fox_was_quickly_being_quick$";
    println!("{}", s.len());
    let a = naive_table(s);
    const B: usize = 8;
    let mut sabtree = SABTree::new(a.len(), B);
    sabtree.build(&a, 0);

    println!("a = {:?}", a);
    println!("btree =\n{:?}", &sabtree.btree);

    let p = "quick}";
    let p0 = "quick";

    let lower_bound = sabtree.search(p0, s);
    let r2 = lower_bound.r;
    println!("{}", &s[r2 as usize..]);

    let upper_bound = sabtree.search(p, s);
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

    let bit = SABTreeIterator::new(&lower_bound, &upper_bound, &sabtree);

    for x in bit {
        println!(
            "matched query {} at {:?}; suffix = {}",
            p0,
            x,
            &s[x.r as usize..]
        );
    }
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

    /*
    let mut ptr = BTreePointer {
        k: 6,
        i: 0,
        r: sabtree.btree[[6, 0]],
    };

    let mut stack = Vec::<BTreePointer>::new();
    println!("{:?}", ptr);
    while ptr.advance(&mut stack, &sabtree) {
        println!("{:?}", ptr);
    }
    */

    Ok(())
}
