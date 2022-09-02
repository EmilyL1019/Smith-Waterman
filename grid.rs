// Adapted from Needleman_Wunsch program on GitHub
struct ImportantExcerpt<'a> {
    part: &'a str,
}
// Creates grid from user entered sequences
pub fn create_grid<'a>(mut seq1: &'a mut String, mut seq2: &'a mut String, len1: i32, len2: i32) -> (Vec<i32>, Vec<String>) {
    // Create new vec and set up row and column labels
    let mut score_grid:Vec<i32>= Vec::new();
    let total_cells = (len1 + 1) * (len2 + 1);
    let mut directions:Vec<String> = Vec::new();
    for _ in 0..total_cells {
        score_grid.push(0);
    }
    // Middle sections
    let mut i = 1;
    while i < total_cells{
        directions = max_cell_score(&mut seq1, &mut seq2, &mut score_grid, &i, &mut directions).to_vec();
        i = i + 1;
    }
    return (score_grid, directions);
}

// Returns direction of the best score for a given cell
fn max_cell_score<'a>(seq1: &'a mut String, seq2: &'a mut String, score_grid: &'a mut Vec<i32>, cell: &'a i32, directions: &'a mut Vec<String>) -> &'a mut Vec<String> {
    // Assign every cell but the first in the first row with L
    if cell > &0 && cell <= &((seq2.len() - 1) as i32) {
        directions.push("L".to_string());
    }
    else if cell % &(seq2.len() as i32) == 0 {
        directions.push("U".to_string());
    }
    else {
        // determine if location is a match
        let seq_1_char_index = ((cell - (cell % seq2.len() as i32)) / seq2.len() as i32 - 1) as usize;
        let seq_2_char_index = (cell % seq2.len() as i32 - 1) as usize;
        let seq_1_char = seq1.chars().nth(seq_1_char_index).unwrap();
        let seq_2_char = seq2.chars().nth(seq_2_char_index).unwrap();
        let mut match_point = 0;
        // Get surrounding scores
        if seq_1_char == seq_2_char {
            // Match scoore
            match_point = 1;
        }
        else {
            match_point = -1            
        }
        // Prevent negative scores
        // Gap score
        let mut from_above:i32 = score_grid[(cell - seq2.len() as i32) as usize] - 1;
        if from_above < 0 {
            from_above = 0;
        }
        let mut from_left:i32 = score_grid[(cell - 1) as usize] - 1;
        if from_left < 0 {
            from_left = 0;
        }
        // Base score
        let mut from_diagonal:i32 = score_grid[(cell - seq2.len() as i32 - 1) as usize] + match_point;
        if from_diagonal < 0 {
            from_diagonal = 0;
        }
        // Find best score
        // Save best directions for each cell
        let mut current_direction = String::new();
        if from_above >= from_left {
            if from_above >= from_diagonal {
                score_grid[*cell as usize] = from_above; 
                current_direction += "U";
            }
            if from_above <= from_diagonal {
                score_grid[*cell as usize] = from_diagonal;
                current_direction += "D";
            }
        }
        if from_above <= from_left {
            if from_left >= from_diagonal{
                score_grid[*cell as usize] = from_left; 
                current_direction += "L";        
            }
            if from_left <= from_diagonal {
                score_grid[*cell as usize] = from_diagonal;
                if !(current_direction.contains("D")) {
                    current_direction += "D";
                }
            }
        }
        directions.push(current_direction);
    }
    return directions;
}