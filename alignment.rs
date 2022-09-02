use std::result;

// Adapted from Needleman_Wunsch program on GitHub
// Transforms the vector of characters into one sequence string
fn char_to_string(characters: Vec<char>) -> String {
    let sequence: String = characters.into_iter().collect();
    return sequence;
}

// Finds all of the cells with the highest score
pub fn highest_cell<'a>(score_grid: &'a Vec<i32>) -> Vec<i32> {
    let mut high_score: i32= 0;
    let mut highest_cells: Vec<i32> = vec![];
    for i in 0..(score_grid.len() - 1) {
        // If higher score is found, empty highest cells vector and put current cell inside
        if score_grid[i] > high_score {
            high_score = score_grid[i];
            highest_cells = vec![i as i32];
        }
        // If another cell with the highest score is found, add to the vector
        else if score_grid[i] == high_score {
            highest_cells.append(&mut vec![i as i32]);            
        }
    }
    // If highest cell score = 0, return -1
    if high_score == 0 {
        highest_cells = vec![-1];
    }
    return highest_cells;
}

// Find best alignment
fn priv_best_alignment<'a>(score_grid: &'a Vec<i32>,directions: &'a Vec<String>, cell: i32, seq1: &'a mut String, seq2: &'a mut String, aligned_seq1: &'a mut Vec<Vec<char>>, aligned_seq2: &'a mut Vec<Vec<char>>, aligned_seq_index: &'a mut usize) -> (Vec<String>, Vec<String>){
    let seq1_char_index = (((cell - (cell % seq2.len() as i32)) / seq2.len() as i32) - 1) as usize;
    let seq2_char_index = ((cell - 1) % seq2.len() as i32) as usize;
    let seq1_char:char = seq1.chars().nth(seq1_char_index).unwrap();
    let seq2_char:char = seq2.chars().nth(seq2_char_index).unwrap();
    let direction_index = (cell - 1) as usize;
    // Recursive cases: Go through directions;
    let mut new_aligned_seq1:&'a mut Vec<Vec<char>> = aligned_seq1;
    let mut new_aligned_seq2: &'a mut Vec<Vec<char>> = aligned_seq2;
    // Keeps track of how many directions the cell has
    let mut has_diagonal:bool = false;
    let mut has_left:bool = false;
    // Look for directions if cell is not 0
    if score_grid[cell as usize] != 0 {
        // If the current direction index is D, add the two corresponding characters to the sequence strings    
        if directions[direction_index].contains("D") {
            let copy_aligned_seq1:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
            let copy_aligned_seq2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
            has_diagonal = true;
            copy_aligned_seq1[*aligned_seq_index].insert(0, seq1_char);
            copy_aligned_seq2[*aligned_seq_index].insert(0, seq2_char);
            // Move to the cell to the diagonally left of the current cell
            if (direction_index as i32 - seq2.len() as i32 - 1) as i32 >= 0 && seq1_char_index > 0 && seq2_char_index > 0{
                priv_best_alignment(score_grid, directions, cell - seq2.len() as i32 - 1, seq1, seq2, copy_aligned_seq1, copy_aligned_seq2, aligned_seq_index);
            }
        }
        // If the current direction index is L, add a gap to the sequence 1 string and the corresponding character to the sequence 2 string
        if directions[direction_index].contains("L") {
            let copy_aligned_seq1:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
            let copy_aligned_seq2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
            has_left = true;
            // Check if cell has multiple paths
            if has_diagonal {
            let new_index = *aligned_seq_index + 1;
                while new_index > *aligned_seq_index {    
                    copy_aligned_seq1.insert(*aligned_seq_index + 1, vec![]);
                    copy_aligned_seq2.insert(*aligned_seq_index + 1,vec![]);
                    *aligned_seq_index = *aligned_seq_index + 1;
                }
            }
            copy_aligned_seq1[*aligned_seq_index].insert(0, '-');
            copy_aligned_seq2[*aligned_seq_index].insert(0, seq2_char);
            // Move to the cell to the left of the current cell
            if (direction_index as i32 - 1) as i32 >= 0 {
                priv_best_alignment(score_grid, directions, cell - 1 as i32, seq1, seq2, copy_aligned_seq1, copy_aligned_seq2, aligned_seq_index);
            }
        }
        // If the current direction index is U, add the corresponding character to the sequence 1 string and a gap to the sequence 2 string
        if directions[direction_index].contains("U") {
            let copy_aligned_seq1_2:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
            let copy_aligned_seq2_2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
            // Check if cell has multiple paths
            if has_diagonal || has_left {
                let new_index = *aligned_seq_index + 1;
                while new_index > *aligned_seq_index {    
                    copy_aligned_seq1_2.insert(*aligned_seq_index + 1, vec![]);
                    copy_aligned_seq2_2.insert(*aligned_seq_index + 1,vec![]);
                    *aligned_seq_index = *aligned_seq_index + 1;
                }
            }
            copy_aligned_seq1_2[*aligned_seq_index].insert(0, seq1_char);
            copy_aligned_seq2_2[*aligned_seq_index].insert(0, '-');    
            // Move to the cell above the current cell
            if (direction_index - seq2.len()) as i32 >= 0 {
                priv_best_alignment(score_grid, directions, cell - seq2.len() as i32, seq1, seq2, copy_aligned_seq1_2, copy_aligned_seq2_2, aligned_seq_index);
            }
        }
    }
    // Base case: Turn finished aligned sequences into strings and return them
    let mut str_aligned_seq1: Vec<String> = Vec::new();
    let mut str_aligned_seq2: Vec<String> = Vec::new();
    for i in 0..(new_aligned_seq1.len() as i32) {
        let char_aligned_seq1:Vec<char> = new_aligned_seq1[i as usize].to_vec().into_iter().collect();
        let char_aligned_seq2:Vec<char> = new_aligned_seq2[i as usize].to_vec().into_iter().collect();
        str_aligned_seq1.append(&mut vec![char_to_string(char_aligned_seq1)]);
        str_aligned_seq2.append(&mut vec![char_to_string(char_aligned_seq2)]);
    }
    return (str_aligned_seq1, str_aligned_seq2);
}
    
fn best_alignment<'a>(score_grid: &'a Vec<i32>,directions: &'a Vec<String>, cell: Vec<i32>, seq1: &'a mut String, seq2: &'a mut String) -> Option<(Vec<String>, Vec<String>)>{
     // if highest cell is -1, return null

   if cell == vec![-1] {
       return None;
    } 
    let mut index:usize = 0;
    let mut str_aligned_seq1:Vec<String> = vec![];
    let mut str_aligned_seq2:Vec<String> = vec![];
    for i in cell {
        let mut aligned_sequence1:Vec<Vec<char>> = vec![[].to_vec()];
        let mut aligned_sequence2:Vec<Vec<char>> = vec![[].to_vec()];
        let (mut group_aligned_seq1, mut group_aligned_seq2) = priv_best_alignment(score_grid, directions, i, seq1, seq2, &mut aligned_sequence1, &mut aligned_sequence2, &mut index);
        str_aligned_seq1.append(&mut group_aligned_seq1);
        str_aligned_seq2.append(&mut group_aligned_seq2);
    }
    for i in 0..str_aligned_seq1.len() {
        let mut j = str_aligned_seq1[i].len();
        while j < str_aligned_seq1[0].len() {
            let char1 = str_aligned_seq1[i - 1].chars().nth(j).unwrap();
            let char2 = str_aligned_seq2[i - 1].chars().nth(j).unwrap();
            str_aligned_seq1[i].insert(j, char1);
            str_aligned_seq2[i].insert(j, char2);
            j = j + 1;
        }
    }
    return Some((str_aligned_seq1, str_aligned_seq2));
}

pub fn build_best_alignment <'a> (score_grid: &'a Vec<i32>,directions: &'a Vec<String>, cell: Vec<i32>, seq1: &'a mut String, seq2: &'a mut String) -> (Vec<String>, Vec<String>){
    let result = best_alignment(score_grid, directions, cell, seq1, seq2);
    match result {
        None => (vec![], vec![]),
        Some(alignments) => alignments,
        
    }
}

// Calculates the score
pub fn score<'a>(str_aligned_seq1: &'a Vec<String>, str_aligned_seq2: &'a Vec<String>) -> i32{
    let mut score: i32 = 0;
    if str_aligned_seq1.len() > 0 {
        for i in 0..str_aligned_seq1[0].len(){
            // Check for match and add 1 to score
            if str_aligned_seq1[0].chars().nth(i) == str_aligned_seq2[0].chars().nth(i) {
                score = score + 1;
            }
            else {
                score = score - 1
            }
        }
    }
    return score;
}

// Prints all alignments
pub fn print_alignments<'a>(str_aligned_seq1: &'a Vec<String>, str_aligned_seq2: &'a Vec<String>, score: i32) {
    if str_aligned_seq1.len() == 0 {
        println!("There are no optimal alignments!");
    }
    else {
        println!("Best Alignments:");
        for i in 0..str_aligned_seq1.len() {
            println!("{}      {}", str_aligned_seq1[i], str_aligned_seq2[i]);
        }
    }
    println!("Score: {}", score);
}
