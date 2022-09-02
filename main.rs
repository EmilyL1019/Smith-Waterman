// Adapted from Needleman_Wunsch program on GitHub
use std::io;
mod grid;
mod alignment;
fn main() {
    // Get two sequences as input
    let mut seq1 : String = String::new();
    io::stdin().read_line(&mut seq1).expect("Cannot read sequence 1.");
    let mut seq2 : String = String::new();
    io::stdin().read_line(&mut seq2).expect("Cannot read sequence 2.");
    let len1 = seq1.len() as i32 - 1;
    let len2 = seq2.len() as i32 - 1;
    let (grid, directions) = grid::create_grid(&mut seq1, &mut seq2, len1, len2);
    let high_cell = alignment::highest_cell(&grid);
    let (aligned_seq1, aligned_seq2) = alignment::build_best_alignment(&grid, &directions, high_cell, &mut seq1, &mut seq2);
    let score = alignment::score(&aligned_seq1, &aligned_seq2);
    alignment::print_alignments(&aligned_seq1, &aligned_seq2, score);
        
}

#[cfg(test)]
mod tests {
    use crate::alignment;
    use crate::alignment::score;
    use crate::grid::create_grid;
    use crate::alignment::build_best_alignment;
    use crate::alignment::print_alignments;
    #[test]
    fn test1() {
        let grid: Vec<i32> = vec![
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0,
        0, 1, 0, 1, 3, 2, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 2, 2, 3, 2, 1, 0, 1,
        0, 0, 0, 0, 1, 1, 3, 4, 3, 2, 1,
        0, 1, 0, 0, 1, 2, 2, 3, 3, 2, 3,
        0, 0, 2, 1, 0, 1, 1, 2, 2, 2, 2,
        0, 0, 1, 3, 2, 1, 0, 1, 3, 3, 2,
        0, 0, 1, 2, 2, 1, 0, 0, 2, 2, 2
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), 
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), 
        "U".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "L".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "D".to_string(), "LD".to_string(), "L".to_string(), "UDL".to_string(), "U".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "UD".to_string(), "D".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "U".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(), "U".to_string(), "UD".to_string(), "U".to_string(), "UD".to_string(), "D".to_string(), "U".to_string(),
        "U".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "D".to_string(), "L".to_string(), 
        "U".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "UD".to_string(), "D".to_string()
        ];
        let mut seq1 : String = "GTCAGGATCT\n".to_string();
        let mut seq2 : String = "ATCAAGGCCA\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 10, 10);
        // Check values
        for i in 0..120 {
            assert_eq!(grid[i], ftn_grid[i]);
        }
        // Check directions
        for i in 0..119 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![73]);
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }
    
    #[test]
    fn test2() {
        let grid = vec![0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 2, 1, 0,
        0, 1, 0, 1, 1, 0,
        0, 0, 0, 0, 2, 2,
        0, 0, 0, 1, 1, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 2, 1,
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(),
        "U".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), "UD".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "LD".to_string()];
        let mut seq1 : String = "ATGCAGGA\n".to_string();
        let mut seq2 : String = "CTGAA\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 8, 5);
        // Check values
        for i in 0..53 {
            assert_eq!(grid[i], ftn_grid[i]);
        }
        // Check directions
        for i in 0..52 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![21, 34, 35, 52]);
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1,&aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["TG", "TGCA", "TGCA", "GA"]);
        assert_eq!(aligned_seq2, vec!["TG", "TG-A", "TGAA", "GA"]);
        assert_eq!(score, 2);
    }
   #[test]
    fn test3() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 2, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 3, 2, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 3, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 1, 1, 1, 2, 2, 1,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 3, 2, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 4, 3, 2, 2, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 3, 3, 2, 3, 2,
        0, 0, 0, 1, 1, 0, 0, 0, 2, 2, 1, 2, 2, 4, 3, 4,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 1, 3, 5, 4,
        0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 3, 2, 4, 4,
        0, 0, 2, 1, 0, 0, 2, 1, 0, 1, 1, 0, 2, 2, 3, 3,
        0, 0, 1, 1, 0, 0, 1, 3, 2, 1, 0, 2, 1, 1, 3, 2,
        0, 0, 1, 0, 0, 0, 1, 2, 2, 3, 2, 1, 1, 0, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 1, 1, 3, 4, 3, 2, 1, 1, 1,
        0, 0, 0, 2, 1, 0, 0, 0, 0, 2, 3, 3, 2, 3, 2, 2,
        0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 2, 4, 3, 2, 4, 3,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 2, 3, 3, 2, 3, 3,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 3, 2, 2, 2, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 2, 2, 1, 1, 1, 1,
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),  "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "LD".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "D".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "D".to_string(), "L".to_string(), "D".to_string(), 
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "LD".to_string(), "U".to_string(), "UD".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UD".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "U".to_string(), "UD".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "LD".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "D".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UD".to_string(), "D".to_string(), "L".to_string(), "UL".to_string(), "UD".to_string(), "D".to_string(), "LD".to_string(), "D".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "U".to_string(), "UD".to_string(), "D".to_string(), "L".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UD".to_string(), "UD".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "L".to_string(),
        "U".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UL".to_string(), "U".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "U".to_string(), "UD".to_string(),
        "U".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "L".to_string(), "UDL".to_string(), "D".to_string(), "UL".to_string(), "UD".to_string(), "D".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), "D".to_string(), "D".to_string(), "LD".to_string(), "UL".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), "UD".to_string(), "D".to_string(), "D".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "U".to_string(), "UD".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "LD".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "D".to_string(), "L".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "U".to_string(), "U".to_string(), "D".to_string(), "L".to_string(), "UL".to_string(), "D".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "U".to_string(), "D".to_string(), "LD".to_string(), "U".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UD".to_string(), "D".to_string(), "UL".to_string(), "UD".to_string(), "D".to_string(), "U".to_string(), "UD".to_string(),
        "U".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UD".to_string(), "D".to_string(), "UDL".to_string(), "UD".to_string(), "UD".to_string(), "UD".to_string()
        ];
        let mut seq1 : String = "AAGTAAGGTGCAGAATGAAA\n".to_string();
        let mut seq2 : String = "CATTCAGGAAGCTGT\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 20, 15);
        //  Check values
        for i in 0..335 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..334 {
            assert_eq!(directions[i], ftn_directions[i]);
        }   
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![174]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AGTAAGGTG"]);
        assert_eq!(aligned_seq2, vec!["AGGAAGCTG"]);
        assert_eq!(score, 5);
    }

    #[test]
    fn test4() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 1, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0,
        0, 0, 0, 0, 0, 1, 0, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "U".to_string(), "D".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        ];
        let mut seq1 : String = "TGACTG\n".to_string();
        let mut seq2 : String = "AAGGTACAA\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 6, 9);
        //  Check values
        for i in 0..69 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..68 {
            assert_eq!(directions[i], ftn_directions[i]);
        }   
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![47]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AC"]);
        assert_eq!(aligned_seq2, vec!["AC"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn test5() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 1,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 2, 1,
        0, 0, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 2, 1,
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), "D".to_string(), 
        "U".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "U".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "L".to_string()
        ];
        let mut seq1 : String = "CTAGATGAG\n".to_string();
        let mut seq2 : String = "TTCAGT\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 9, 6);
        //  Check values
        for i in 0..69 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..68 {
            assert_eq!(directions[i], ftn_directions[i]);
        }   
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![33, 48, 68]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AG", "AGAT", "AG"]);
        assert_eq!(aligned_seq2, vec!["AG", "AG-T", "AG"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn test6() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "D".to_string(), "D".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "D".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "D".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        ];
        let mut seq1 : String = "TTGATGT\n".to_string();
        let mut seq2 : String = "AAACTACA\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 8, 8);
        //  Check values
        for i in 0..63{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..62 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![14, 23, 37, 38, 39, 42, 44, 50, 68]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["T", "T", "A", "A", "A", 
        "A", "A", "T", "T"]);
        assert_eq!(aligned_seq2, vec!["T", "T", "A", "A", "A", 
        "A", "A", "T", "T"]);
        assert_eq!(score, 1);
    }

    #[test]
    fn test7() {
        let grid = vec![0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0
        ];
        let directions:Vec<String> = vec!["L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), "L".to_string(), 
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(),
        "U".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string(), "UDL".to_string()
        ];
        let mut seq1 : String = "TTTTT\n".to_string();
        let mut seq2 : String = "AAAAA\n".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 5, 5);
        //  Check values
        for i in 0..35{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..34 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = alignment::highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![-1]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(score, 0);
    }
}
