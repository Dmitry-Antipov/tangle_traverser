#!/usr/bin/env python3

import sys
import difflib
#TODO:LLM-generated, way too much excessive

def highlight_string_differences(str1, str2, max_width=80):
    """
    Highlight character-level differences between two strings.
    Shows exact positions where strings differ.
    """
    if str1 == str2:
        return "Strings are identical"
    
    result = []
    
    # Find all differing positions
    min_len = min(len(str1), len(str2))
    max_len = max(len(str1), len(str2))
    diff_positions = []
    
    # Check character by character
    for i in range(min_len):
        if str1[i] != str2[i]:
            diff_positions.append(i)
    
    # Add positions for length differences
    if len(str1) != len(str2):
        for i in range(min_len, max_len):
            diff_positions.append(i)
    
    if not diff_positions:
        return "Strings are identical"
    
    result.append(f"String lengths: {len(str1)} vs {len(str2)}")
    result.append(f"Differences found at {len(diff_positions)} position(s): {diff_positions[:10]}{'...' if len(diff_positions) > 10 else ''}")
    
    # Show detailed diff for each differing segment
    for i, pos in enumerate(diff_positions[:5]):  # Show first 5 differences
        start = max(0, pos - 10)
        end = min(max_len, pos + 11)
        
        result.append(f"\nDifference {i+1} at position {pos}:")
        
        # Show context around the difference
        context1 = str1[start:end] if pos < len(str1) else str1[start:] + "<END>"
        context2 = str2[start:end] if pos < len(str2) else str2[start:] + "<END>"
        
        result.append(f"  String 1 [{start}:{end}]: {repr(context1)}")
        result.append(f"  String 2 [{start}:{end}]: {repr(context2)}")
        
        # Create visual pointer to the exact difference
        if pos < len(str1) and pos < len(str2):
            pointer_pos = pos - start
            if pointer_pos >= 0 and pointer_pos < len(context1):
                pointer1 = ' ' * (19 + pointer_pos) + '^'
                pointer2 = ' ' * (19 + pointer_pos) + '^'
                result.append(f"  {pointer1} (char: {repr(str1[pos]) if pos < len(str1) else 'EOF'})")
                result.append(f"  {pointer2} (char: {repr(str2[pos]) if pos < len(str2) else 'EOF'})")
        elif pos >= len(str1):
            result.append(f"  String 1 ends here, String 2 continues with: {repr(str2[pos]) if pos < len(str2) else 'EOF'}")
        elif pos >= len(str2):
            result.append(f"  String 2 ends here, String 1 continues with: {repr(str1[pos]) if pos < len(str1) else 'EOF'}")
    
    if len(diff_positions) > 5:
        result.append(f"\n... and {len(diff_positions) - 5} more differences")
    
    # Add unified diff for overall context
    result.append("\nUnified diff:")
    diff = list(difflib.unified_diff(
        str1.splitlines(keepends=True), 
        str2.splitlines(keepends=True),
        fromfile="String 1",
        tofile="String 2",
        lineterm=""
    ))
    
    if diff:
        result.extend(diff[:20])  # Limit diff output
    
    return '\n'.join(result)

def compare_files(file1, file2):
    """Compare two files line by line and highlight differences."""
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # Use difflib for a comprehensive comparison
    diff = list(difflib.unified_diff(
        lines1, 
        lines2,
        fromfile=file1,
        tofile=file2,
        lineterm=""
    ))
    
    if not diff:
        print("Files are identical")
        return
    
    print("File differences:")
    for line in diff:
        print(line)
    
    # Also check for line-by-line string differences for detailed analysis
    max_lines = max(len(lines1), len(lines2))
    detailed_diffs = False
    
    for i in range(max_lines):
        line1 = lines1[i].rstrip('\n') if i < len(lines1) else ""
        line2 = lines2[i].rstrip('\n') if i < len(lines2) else ""
        
        if line1 != line2:
            if not detailed_diffs:
                print("\n" + "="*50)
                print("DETAILED LINE-BY-LINE ANALYSIS:")
                print("="*50)
                detailed_diffs = True
            
            print(f"\nLine {i+1} differs:")
            diff_result = highlight_string_differences(line1, line2)
            print(diff_result)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python diff_compare.py file1 file2")
        sys.exit(1)
    
    compare_files(sys.argv[1], sys.argv[2])

