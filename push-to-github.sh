#!/bin/bash
# Script to push ViroAgent Analysis repository to GitHub
# Usage: ./push-to-github.sh [github-username] [repository-name]

set -e

# Default values
GITHUB_USER="${1:-your-username}"
REPO_NAME="${2:-viroagent-analysis}"

echo "============================================="
echo "ViroAgent Analysis - GitHub Repository Setup"
echo "============================================="
echo ""

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo "❌ Git is not installed. Please install git first."
    exit 1
fi

echo "📦 Repository: $REPO_NAME"
echo "👤 GitHub User: $GITHUB_USER"
echo ""

# Step 1: Check remote
if git remote | grep -q origin; then
    echo "✅ Remote 'origin' already configured."
    CURRENT_REMOTE=$(git remote get-url origin)
    echo "   Current remote URL: $CURRENT_REMOTE"
    echo ""
else
    echo "🔧 Step 1: Configuring remote repository..."
    echo ""
    echo "Please choose an authentication method:"
    echo "1) SSH (requires SSH keys configured on GitHub)"
    echo "2) HTTPS (requires personal access token)"
    echo ""
    read -p "Enter choice (1 or 2): " AUTH_CHOICE
    
    if [ "$AUTH_CHOICE" = "1" ]; then
        REMOTE_URL="git@github.com:$GITHUB_USER/$REPO_NAME.git"
        echo "Using SSH URL: $REMOTE_URL"
        echo ""
        echo "⚠️  Make sure you have:"
        echo "   - SSH key added to your GitHub account"
        echo "   - Repository '$REPO_NAME' created on GitHub"
    elif [ "$AUTH_CHOICE" = "2" ]; then
        REMOTE_URL="https://github.com/$GITHUB_USER/$REPO_NAME.git"
        echo "Using HTTPS URL: $REMOTE_URL"
        echo ""
        echo "⚠️  Make sure you have:"
        echo "   - Repository '$REPO_NAME' created on GitHub"
        echo "   - GitHub personal access token with 'repo' scope"
        echo "   (You'll be prompted for credentials when pushing)"
    else
        echo "❌ Invalid choice. Exiting."
        exit 1
    fi
    
    echo ""
    read -p "Have you created the repository '$REPO_NAME' on GitHub? (y/n): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo ""
        echo "📝 Please create the repository first:"
        echo "   1. Go to https://github.com/new"
        echo "   2. Repository name: '$REPO_NAME'"
        echo "   3. Description: 'Computational virology research assistant'"
        echo "   4. Choose public or private"
        echo "   5. Click 'Create repository'"
        echo ""
        read -p "Press Enter after creating the repository..."
    fi
    
    git remote add origin "$REMOTE_URL"
    echo "✅ Remote 'origin' added."
fi

echo ""
echo "🔍 Step 2: Checking repository status..."
git status

echo ""
echo "🚀 Step 3: Pushing to GitHub..."
echo "   Branch: main"
echo "   Remote: origin"
echo ""

# Rename master to main if needed
if git branch | grep -q "master"; then
    if ! git branch | grep -q "main"; then
        echo "🔄 Renaming branch from 'master' to 'main'..."
        git branch -m main
    fi
fi

# Push with force if needed (first push)
echo "📤 Pushing code..."
if git push -u origin main 2>&1 | tee /tmp/push_output; then
    echo ""
    echo "✅ Successfully pushed to GitHub!"
    echo ""
    echo "🌐 Your repository is now available at:"
    echo "   https://github.com/$GITHUB_USER/$REPO_NAME"
else
    echo ""
    echo "⚠️  Push may have failed. Common issues:"
    echo "   - Repository not created on GitHub"
    echo "   - Authentication problems (SSH keys or token)"
    echo "   - Network connectivity"
    echo ""
    echo "🔧 Try these steps:"
    echo "   1. Create repository at https://github.com/new"
    echo "   2. Ensure your SSH key is added to GitHub (for SSH)"
    echo "   3. For HTTPS, use: git push https://<token>@github.com/$GITHUB_USER/$REPO_NAME.git"
    echo ""
    exit 1
fi

echo ""
echo "============================================="
echo "🎉 Setup Complete!"
echo "============================================="
echo ""
echo "Next steps:"
echo "1. Review the repository on GitHub"
echo "2. Set up GitHub Pages for documentation (optional)"
echo "3. Share with collaborators"
echo "4. Star the repository if you find it useful!"
echo ""
echo "For issues or questions, check the README.md file."